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
		virtual unsigned int requiredADdirs() const CADET_NOEXCEPT; // ?

		virtual UnitOpIdx unitOperationId() const CADET_NOEXCEPT { return _unitOpIdx; }
		virtual unsigned int numComponents() const CADET_NOEXCEPT { return _disc.nComp; }
		virtual void setFlowRates(active const* in, active const* out) CADET_NOEXCEPT; // ?
		virtual unsigned int numInletPorts() const CADET_NOEXCEPT { return 1; }
		virtual unsigned int numOutletPorts() const CADET_NOEXCEPT { return 1; }
		virtual bool canAccumulate() const CADET_NOEXCEPT { return false; } // ?

		static const char* identifier() { return "LUMPED_RATE_MODEL_WITHOUT_PORES_DG"; }
		virtual const char* unitOperationName() const CADET_NOEXCEPT { return identifier(); }

		virtual bool configureModelDiscretization(IParameterProvider& paramProvider, IConfigHelper& helper);
		virtual bool configure(IParameterProvider& paramProvider);
		virtual void notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const ConstSimulationState& simState, const AdJacobianParams& adJac);

		virtual void useAnalyticJacobian(const bool analyticJac);

		virtual void reportSolution(ISolutionRecorder& recorder, double const* const solution) const;
		virtual void reportSolutionStructure(ISolutionRecorder& recorder) const; // modify ?

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
		// lean ?
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
		// modify ?
		virtual void multiplyWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* yS, double alpha, double beta, double* ret);
		virtual void multiplyWithDerivativeJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* sDot, double* ret);

		inline void multiplyWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* yS, double* ret)
		{
			multiplyWithJacobian(simTime, simState, yS, 1.0, 0.0, ret);
		}

		virtual bool setParameter(const ParameterId& pId, double value);
		virtual bool setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue);
		virtual void setSensitiveParameterValue(const ParameterId& id, double value);
		// modify ?
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
		// modify !
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
		{ 	// NOTE: no different Riemann solvers, BC, BF compared to external implementation
		public:
			unsigned int nComp; //!< Number of components
			unsigned int nCol; //!< Number of column cells
			unsigned int polyDeg; //!< polynomial degree
			unsigned int nNodes; //!< Number of nodes per cell
			unsigned int nPoints; //!< Number of discrete Points
			bool modal;	//!< bool switch: 1 for modal basis, 0 for nodal basis
			Eigen::VectorXd nodes; //!< Array with positions of nodes in reference element
			Eigen::MatrixXd polyDerM; //!< Array with polynomial derivative Matrix
			Eigen::VectorXd invWeights; //!< Array with weights for numerical quadrature of size nNodes
			Eigen::MatrixXd invMM; //!< dense !INVERSE! mass matrix for modal (exact) integration

			unsigned int* nBound; //!< Array with number of bound states for each component
			unsigned int* boundOffset; //!< Array with offset to the first bound state of each component in the solid phase
			unsigned int strideBound; //!< Total number of bound states

			// vgl. convDispOperator
			std::vector<double> dispersion; //!< Column dispersion (may be section dependent) \f$ D_{\text{ax}} \f$
			std::vector<double> velocity; //!< Interstitial velocity (may be section dependent) \f$ u \f$
			double crossSection; //!< Cross section area 
			double length_;
			double deltaZ;
			double porosity;

			//
			std::vector<Eigen::MatrixXd> coefs;
			unsigned int nSec;

			//@TODO ConvDisp  DG operator, cache
			Eigen::VectorXd g;
			Eigen::VectorXd h;
			Eigen::VectorXd surfaceFlux; //!< stores the surface flux values
			Eigen::Vector4d boundary; //!< stores the boundary values from Danckwert boundary conditions

			//@TODO: unnecessary when binding model implementation is used
			std::vector<double> adsorption;
			std::vector<double> desorption;
			std::vector<double> ADratio;
			std::vector<bool> isKinetic;
			std::string isotherm;

			//Eigen::SparseMatrix<double> _JacC; //!< Static Convection Dispersion Jacobian
			Eigen::MatrixXd _JacC;

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

				//@TODO ConvDisp DG Operator
				g.resize(nPoints);
				g.setZero();
				h.resize(nPoints);
				h.setZero();
				boundary.setZero();
				surfaceFlux.resize(nCol + 1);
				surfaceFlux.setZero();

				// @TODO: binding model!
				isKinetic.resize(nComp);
				for (int comp = 0; comp < nComp; comp++) {
					//coefs[comp] = MatrixXd::Zero(nSec, 4);
					isKinetic[comp] = false;
				}

				lglNodesWeights();
				// @TODO: make modal/nodal switch during calculation possible?
				if (modal) {
					derivativeMatrix_JACOBI();
					invMMatrix_JACOBI();
				}
				else {
					// @TODO: check why hdf5 cant be closed when D matrix is determined with lagrange..
					derivativeMatrix_LAGRANGE();
					//derivativeMatrix_JACOBI(); // D matrix is approximately the same
					invMMatrix_JACOBI();
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
					for (int j = 1; j <= floor((polyDeg + 1) / 2) - 1; j++) {
						//  first guess for Newton iteration
						nodes[j] = -cos(M_PI * (j + 0.25) / polyDeg - 3 / (8.0 * polyDeg * M_PI * (j + 0.25)));
						// Newton iteration to find roots of Legendre Polynomial
						for (int k = 0; k <= nIterations; k++) {
							qAndL(nodes[j], L, q, qder);
							nodes[j] = nodes[j] - q / qder;
							if (abs(q / qder) <= tolerance * abs(nodes[j])) {
								break;
							}
						}
						// calculate weights
						qAndL(nodes[j], L, q, qder);
						invWeights[j] = 2 / (polyDeg * (polyDeg + 1.0) * pow(L, 2));
						nodes[polyDeg - j] = -nodes[j]; // copy to second half of points and weights
						invWeights[polyDeg - j] = invWeights[j];
					}
				}
				if (polyDeg % 2 == 0) { // for even polyDeg we have an odd number of points which include 0.0
					qAndL(0.0, L, q, qder);
					nodes[polyDeg / 2] = 0;
					invWeights[polyDeg / 2] = 2 / (polyDeg * (polyDeg + 1.0) * pow(L, 2));
				}
				// inverse the weights
				invWeights = invWeights.cwiseInverse();
			}

			/**
			 * @brief computation of barycentric weights for fast polynomial evaluation
			 * @param [in] weights pre-allocated vector for barycentric weights. Must be set to ones!
			 */
			void barycentricWeights(VectorXd& baryWeights) {
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
			 * @brief computation of LAGRANGE polynomial derivative matrix
			 */
			void derivativeMatrix_LAGRANGE() {
				VectorXd baryWeights = VectorXd::Ones(polyDeg + 1);
				barycentricWeights(baryWeights);
				for (int i = 0; i <= polyDeg; i++) {
					for (int j = 0; j <= polyDeg; j++) {
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
				for (int i = 0; i < nodes.size(); ++i) {
					PL(0, i) = 1.0 / std::sqrt(gamma0);
				}
				if (N == 0) {
					for (int i = 0; i < nodes.size(); ++i) {
						P(i, index) = PL(0, i);
					}
					return;
				}
				double gamma1 = (alpha + 1) * (beta + 1) / (alpha + beta + 3) * gamma0;
				for (int i = 0; i < nodes.size(); ++i) {
					PL(1, i) = ((alpha + beta + 2) * nodes(i) / 2 + (alpha - beta) / 2) / std::sqrt(gamma1);
				}
				if (N == 1) {
					for (int i = 0; i < nodes.size(); ++i) {
						P(i, index) = PL(1, i);
					}
					return;
				}
				double a = 2.0 / (2.0 + alpha + beta) * std::sqrt((alpha + 1) * (beta + 1) / (alpha + beta + 3));
				for (int i = 0; i < N - 1; ++i) {
					double j = i + 1.0;
					double h1 = 2.0 * (i + 1.0) + alpha + beta;
					double a_aux = 2.0 / (h1 + 2.0) * std::sqrt((j + 1.0) * (j + 1.0 + alpha + beta) * (j + 1.0 + alpha) * (j + 1.0 + beta) / (h1 + 1) / (h1 + 3));
					double b = -(std::pow(alpha, 2) - std::pow(beta, 2)) / h1 / (h1 + 2);

					for (int k = 0; k < nodes.size(); ++k) {
						PL(i + 2, k) = 1 / a_aux * (-a * PL(i, k) + (nodes(k) - b) * PL(i + 1, k));
					}
					a = a_aux;
				}
				for (int i = 0; i < nodes.size(); ++i) {
					P(i, index) = PL(N, i);
				}
			}

			/**
			* @brief returns generalized Jacobi Vandermonde matrix
			* @detail normalized Legendre Vandermonde matrix
			*/
			MatrixXd getVandermonde_JACOBI() {

				MatrixXd V(nodes.size(), nodes.size());

				for (int j = 0; j < V.cols(); ++j) {
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

				for (int order = 1; order < nodes.size(); order++) {
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
				for (int i = 0; i < Np; ++i) {
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

	//	IExternalFunction* _extFun; //!< External function (owned by library user) @TODO: what for, inlet?

		// TODO: DG ConvDisp implementation
		parts::ConvectionDispersionOperatorBase _convDispOp; //!< Convection dispersion operator for interstitial volume transport

		Eigen::MatrixXd _jac; //!< Jacobian
		Eigen::MatrixXd _jacDisc; //!< Jacobian with time derivatives from BDF method
		Eigen::MatrixXd _jacInlet; //!< Jacobian inlet DOF block matrix connects inlet DOFs to first bulk cells

		active _totalPorosity; //!< Total porosity \f$ \varepsilon_t \f$

		bool _analyticJac; //!< Determines whether AD or analytic Jacobians are used
		unsigned int _jacobianAdDirs; //!< Number of AD seed vectors required for Jacobian computation
		// TODO: factorization for DG jacobian?
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
			inline int strideColCell() const CADET_NOEXCEPT { return static_cast<int>(_disc.nNodes); }
			inline int strideColNode() const CADET_NOEXCEPT { return 1; }
			inline int strideColComp() const CADET_NOEXCEPT { return static_cast<int>(_disc.nNodes * _disc.nCol); }

			// @TODO: modify !
			inline int strideColLiquid() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp); }
			inline int strideColBound() const CADET_NOEXCEPT { return static_cast<int>(_disc.strideBound); }
			// @TODO: modify !
			// Offsets
			inline int offsetC() const CADET_NOEXCEPT { return _disc.nComp; }
			inline int offsetBoundComp(unsigned int comp) const CADET_NOEXCEPT { return _disc.boundOffset[comp]; }
			// @TODO: modify!
			// Return pointer to first element of state variable in state vector
			template <typename real_t> inline real_t* c(real_t* const data) const { return data + offsetC(); }
			template <typename real_t> inline real_t const* c(real_t const* const data) const { return data + offsetC(); }

			template <typename real_t> inline real_t* q(real_t* const data) const { return data + offsetC() + strideColLiquid(); }
			template <typename real_t> inline real_t const* q(real_t const* const data) const { return data + offsetC() + strideColLiquid(); }
			// @TODO: modify!
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
			// @TODO?? actually we need number of axial discrete points here, not number of axial cells !
			virtual unsigned int numAxialCells() const CADET_NOEXCEPT { return _disc.nCol * _disc.nNodes; }
			virtual unsigned int numRadialCells() const CADET_NOEXCEPT { return 0; }
			virtual unsigned int numInletPorts() const CADET_NOEXCEPT { return 1; }
			virtual unsigned int numOutletPorts() const CADET_NOEXCEPT { return 1; }
			virtual unsigned int numParticleTypes() const CADET_NOEXCEPT { return 1u; }
			virtual unsigned int numParticleShells(unsigned int parType) const CADET_NOEXCEPT { return 0u; }
			virtual unsigned int numBoundStates(unsigned int parType) const CADET_NOEXCEPT { return _disc.strideBound; }
			virtual unsigned int numBulkDofs() const CADET_NOEXCEPT { return _disc.nComp * _disc.nCol * _disc.nNodes; }
			virtual unsigned int numParticleMobilePhaseDofs(unsigned int parType) const CADET_NOEXCEPT { return 0u; }
			virtual unsigned int numSolidPhaseDofs(unsigned int parType) const CADET_NOEXCEPT { return _disc.strideBound * _disc.nCol * _disc.nNodes; }
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

			/**
			*@brief calculates the physical node coordinates of DG discretization with double! interface nodes
			*/
			virtual void axialCoordinates(double* coords) const {
				VectorXd x_l = VectorXd::LinSpaced(_disc.nCol + 1u, 0.0, _disc.length_);
				for (int i = 0; i < _disc.nCol; i++) {
					for (int j = 0; j < _disc.nNodes; j++) {
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



		// ==========================================================================================================================================================  //
		// ==========================================================================================================================================================  //
		// ========================================						DG Residual, RHS						======================================================  //
		// ==========================================================================================================================================================  //
		// ==========================================================================================================================================================  //

		/**
* @brief calculates the volume Integral of the auxiliary equation
* @param [in] current state vector
* @param [in] stateDer vector to be changed
* @param [in] aux true if auxiliary, else main equation
*/
		void volumeIntegral(Discretization& disc, const double* const statePtr, double* const stateDerPtr) {
			// Map pointer to Eigen Objects
			Eigen::Map<const Eigen::VectorXd> state(statePtr, disc.nNodes * disc.nCol);
			Eigen::Map<Eigen::VectorXd> stateDer(stateDerPtr, disc.nNodes * disc.nCol);

			// comp-cell-node state vector: use of Eigen lib performance
			for (int Cell = 0; Cell < disc.nCol; Cell++) {
				stateDer.segment(Cell * disc.nNodes, disc.nNodes)
					-= disc.polyDerM * state.segment(Cell * disc.nNodes, disc.nNodes);
			}
		}

		/*
		* @brief calculates the interface fluxes h* of Convection Dispersion equation
		*/
		void InterfaceFlux(const double* const cPtr, const double* const gPtr, Discretization& disc, unsigned int secIdx) {
			// Map pointer to Eigen Objects
			Eigen::Map<const Eigen::VectorXd> C(cPtr, disc.nNodes * disc.nCol);
			Eigen::Map<const Eigen::VectorXd> g(gPtr, disc.nNodes * disc.nCol);

			Indexer idxr(disc);
			int strideCell = idxr.strideColCell();
			int strideNode = idxr.strideColNode();

			// Conv.Disp. flux: h* = h*_conv + h*_disp = numFlux(v c_l, v c_r) + 0.5 sqrt(D_ax) (S_l + S_r)
			
			// calculate inner interface fluxes
			for (int Cell = 1; Cell < disc.nCol; Cell++) {
				// h* = h*_conv + h*_disp
				disc.surfaceFlux[Cell] // inner interfaces
					= disc.velocity[secIdx] * (C[Cell * strideCell - strideNode])
					- 0.5 * std::sqrt(disc.dispersion[secIdx]) * (g[Cell * strideCell - strideNode] // left cell
						                                          + g[Cell * strideCell]);
			}

			// boundary fluxes
				// left boundary interface
			disc.surfaceFlux[0]
				= disc.velocity[secIdx] * disc.boundary[0];

			// right boundary interface
			disc.surfaceFlux[disc.nCol]
				= disc.velocity[secIdx] * (C[disc.nCol * strideCell - strideNode])
				- std::sqrt(disc.dispersion[secIdx]) * 0.5 * (g[disc.nCol * strideCell - strideNode] // last cell last node
					+ disc.boundary[3]); // right boundary value S
		}

		/**
		* @brief calculates and fills the surface flux values for auxiliary equation
		* @param [in] aux true if auxiliary, else main equation
		*/
		void InterfaceFluxAuxiliary(const double* const cPtr, Discretization& disc) {

			// Map pointer to Eigen Objects
			Eigen::Map<const Eigen::VectorXd> C(cPtr, disc.nNodes * disc.nCol);

			Indexer idxr(disc);
			int strideCell = idxr.strideColCell();
			int strideNode = idxr.strideColNode();

			// Auxiliary flux: c* = 0.5 (c_l + c_r)

			// calculate inner interface fluxes
			for (int Cell = 1; Cell < disc.nCol; Cell++) {
				disc.surfaceFlux[Cell] // left interfaces
					= 0.5 * (C[Cell * strideCell - strideNode] + // left node
						C[Cell * strideCell]); // right node
			}
			// calculate boundary interface fluxes

			disc.surfaceFlux[0] // left boundary interface
				= 0.5 * (C[0] + // boundary value
					C[0]); // first cell first node

			disc.surfaceFlux[(disc.nCol)] // right boundary interface
				= 0.5 * (C[disc.nCol * strideCell - strideNode] + // last cell last node
					C[disc.nCol * strideCell - strideNode]);// // boundary value
		}

		/**
		* @brief calculates the surface Integral, depending on the approach (modal/nodal)
		* @param [in] state relevant state vector
		* @param [in] stateDer state derivative vector the solution is added to
		* @param [in] aux true for auxiliary, false for main equation
			surfaceIntegral(cPtr, &(disc.g[0]), disc,&(disc.h[0]), resPtrC, 0, secIdx);
		*/ // 
		void surfaceIntegral(const double* const cPtr, const double* const gPtr, Discretization& disc, const double* const statePtr, double* const stateDerPtr, bool aux, unsigned int secIdx) {

			Eigen::Map<const Eigen::VectorXd> state(statePtr, disc.nNodes * disc.nCol);
			Eigen::Map<Eigen::VectorXd> stateDer(stateDerPtr, disc.nNodes * disc.nCol);

			Indexer idxr(disc);
			int strideCell = idxr.strideColCell();
			int strideNode = idxr.strideColNode();

			// calc numerical flux values c* or h* depending on equation switch aux
			(aux == 1) ? InterfaceFluxAuxiliary(cPtr, disc) : InterfaceFlux(cPtr, &(disc.g[0]), disc, secIdx);
			if (disc.modal) { // modal approach -> dense mass matrix
				for (int Cell = 0; Cell < disc.nCol; Cell++) {
					// strong surface integral -> M^-1 E [state - state*]
					for (int Node = 0; Node < disc.nNodes; Node++) {
						stateDer[Cell * strideCell + Node * strideNode]
							-= disc.invMM(Node, 0) * (state[Cell * strideCell]
								- disc.surfaceFlux[Cell])
							- disc.invMM(Node, disc.polyDeg) * (state[Cell * strideCell + disc.polyDeg * strideNode]
								- disc.surfaceFlux[(Cell + 1)]);
					}
				}
			}
			else { // nodal approach -> diagonal mass matrix
				for (int Cell = 0; Cell < disc.nCol; Cell++) {
					// strong surface integral -> M^-1 B [state - state*]
					stateDer[Cell * strideCell /*+ Comp * strideComp*/] // first cell node
						-= disc.invWeights[0] * (state[Cell * strideCell] // first node
							- disc.surfaceFlux(Cell));
					stateDer[Cell * strideCell + disc.polyDeg * strideNode] // last cell node
						+= disc.invWeights[disc.polyDeg] * (state[Cell * strideCell + disc.polyDeg * strideNode]
							- disc.surfaceFlux(Cell + 1));
				}
			}
		}

		/**
		* @brief calculates the substitute h = vc - sqrt(D_ax) S(c)
		*/
		void calcH(const double* const cPtr, Discretization& disc, unsigned int secIdx) {
			// Map pointer to Eigen Objects
			Eigen::Map<const Eigen::VectorXd> C(cPtr, disc.nNodes * disc.nCol);
			disc.h = disc.velocity[secIdx] * C - std::sqrt(disc.dispersion[secIdx]) * disc.g;
		}

		/*
		* calculates isotherm RHS
		*/
		void calcRHSq_DG(const double* const yPtr, double* const resPtr, Discretization& disc) {

			// Map pointer to Eigen Objects
			Eigen::Map<const Eigen::VectorXd> C(yPtr + disc.nComp, disc.nComp * disc.nNodes * disc.nCol);
			Eigen::Map<const Eigen::VectorXd> q(yPtr + disc.nComp + disc.nComp * disc.nNodes * disc.nCol, disc.nComp * disc.nNodes * disc.nCol);
			Eigen::Map<Eigen::VectorXd> qRes(resPtr + disc.nComp + disc.nComp * disc.nPoints, disc.nComp * disc.nPoints);

			Indexer idxr(disc);
			int strideCell = idxr.strideColCell();
			int strideComp = idxr.strideColComp();
			int strideNode = idxr.strideColNode();

			// reset cache
			qRes.setZero();

			for (int comp = 0; comp < disc.nComp; comp++) {
				if (disc.isotherm == "Linear" && disc.isKinetic[comp]) {
					qRes.segment(comp * strideComp, strideComp) = disc.adsorption[comp] * C.segment(comp * strideComp, strideComp)
						- disc.desorption[comp] * q.segment(comp * strideComp, strideComp);
				}
				else if (disc.isotherm == "Linear" && !disc.isKinetic[comp]) { // NOTE: equilibrium constant stored in adsorption
					qRes.segment(comp * strideComp, strideComp) = disc.adsorption[comp] * C.segment(comp * strideComp, strideComp)
						- q.segment(comp * strideComp, strideComp);
				}
				else if (disc.isotherm == "Langmuir" && !disc.isKinetic[comp]) {
					double factor;
					int nPoints = disc.nCol * disc.nNodes * disc.nComp;
					for (int point = 0; point < nPoints; point++) {
						factor = 1.0;
						for (int nComp = 0; nComp < disc.nComp; nComp++) {
							// calc(1 + sum(b_i * c_i))
							for (int comp = 0; comp < disc.nComp; comp++) {
								factor += disc.ADratio[comp] * C[point * strideNode + comp * strideComp];
							}
							qRes[point * strideNode + nComp * strideComp]
								= (disc.adsorption[nComp] * C[point * strideNode + nComp * strideComp]) / factor
								- q[point * strideNode + nComp * strideComp];
						}
					}
				}
			}
		}

		/**
		* @brief applies the inverse Jacobian of the mapping
		*/
		void applyMapping(const Discretization& disc, double* const statePtr) {
			Eigen::Map<Eigen::VectorXd> state(statePtr, disc.nPoints);
			state = state * (2 / disc.deltaZ);
		}
		/**
		* @brief applies the inverse Jacobian of the mapping and Auxiliary factor
		*/
		void applyMapping_Aux(const Discretization& disc, double* const statePtr, unsigned int secIdx) {
			Eigen::Map<Eigen::VectorXd> state(statePtr, disc.nPoints);
			state = state * (-2 / disc.deltaZ) * ((disc.dispersion[secIdx] == 0) ? 1.0 : std::sqrt(disc.dispersion[secIdx]));
		}
		/**
		* @brief calculates q from c (equilibrium ! )
		* @param nPoints number of discrete points -> state.size()
		*/
		void calcQ(const VectorXd& state, VectorXd& q, const Discretization& disc, unsigned int secIdx) {

			Indexer idxr(disc);
			int strideCell = idxr.strideColCell();
			int strideComp = idxr.strideColComp();
			int strideNode = idxr.strideColNode();

			if (disc.isotherm == "Linear" && !disc.isKinetic[secIdx]) {
				for (int nComp = 0; nComp < disc.nComp; nComp++) {
					q.segment(nComp * strideComp, strideComp) = disc.adsorption[nComp] *
						state.segment(nComp * strideComp, strideComp);
				}
			}
			else if (disc.isotherm == "Langmuir" && !disc.isKinetic[secIdx]) {
				double factor;
				for (int point = 0; point < disc.nPoints; point++) {
					factor = 1.0;
					for (int nComp = 0; nComp < disc.nComp; nComp++) {
						// calc(1 + sum(b_i * c_i))
						for (int comp = 0; comp < disc.nComp; comp++) {
							factor += disc.ADratio[comp] * state[point * strideNode + comp * strideComp];
						}
						q[point * strideNode + nComp * strideComp]
							= (disc.adsorption[nComp] * state[point * strideNode + nComp * strideComp]) / factor;
					}
				}
			}
			else if (disc.isotherm == "Linear" && disc.isKinetic[secIdx]) {
				for (int nComp = 0; nComp < disc.nComp; nComp++) {
					q.segment(nComp * strideComp, strideComp) = disc.adsorption[nComp] * state.segment(nComp * strideComp, strideComp)
														      - disc.desorption[nComp] * q.segment(nComp * strideComp, strideComp);
				}
			}
			else {
				throw std::invalid_argument("spelling error or this isotherm is not implemented yet");
			}
		}

		void ConvDisp_DG(const double* const cPtr, double* const resPtrC, Discretization& disc, double t, unsigned int Comp, unsigned int secIdx) {

			// Map pointer to Eigen Objects
			Eigen::Map<const Eigen::VectorXd> C(cPtr, disc.nPoints);
			Eigen::Map<Eigen::VectorXd> resC(resPtrC, disc.nPoints);

			// ===================================//
			// reset cache                       //
			// =================================//

			resC.setZero();
			disc.h.setZero();
			disc.g.setZero();
			disc.surfaceFlux.setZero();

			// ==================================//
			// solve auxiliary system g = dc/dx  //
			// ==================================//

			// DG volumne and surface integral in strong form
			volumeIntegral(disc, cPtr, &(disc.g[0]));
			// surface integral in strong form
			surfaceIntegral(cPtr, &(disc.g[0]), disc, cPtr, &(disc.g[0]), 1, secIdx);
			// inverse mapping from reference space and auxiliary factor
			applyMapping_Aux(disc, &(disc.g[0]), secIdx);
			// reset surface flux storage as it is used twice
			std::cout << "AUXSurfFlux: \n" << disc.surfaceFlux << std::endl;
			disc.surfaceFlux.setZero();

			// ===================================//
			// solve main equation w_t = dh/dx   //
			// =================================//

			// calculate the substitute h(S(c), c) = sqrt(D_ax) S(c) - v c
			calcH(cPtr, disc, secIdx);
			// DG volumne and surface integral in strong form
			volumeIntegral(disc, &(disc.h[0]), resPtrC);
			// update boundary values ( with S)
			calcBoundaryValues(disc, cPtr, Comp, secIdx);
			surfaceIntegral(cPtr, &(disc.g[0]), disc,&(disc.h[0]), resPtrC, 0, secIdx);
			std::cout << "MAINSurfFlux: \n" << disc.surfaceFlux << std::endl;

			// inverse mapping from reference space
			applyMapping(disc, resPtrC);

		}


		void calcBoundaryValues(Discretization& disc, const double* const c, unsigned int comp, unsigned int secIdx) {

			Indexer idxr(disc);
			int strideCell = idxr.strideColCell();
			int strideComp = idxr.strideColComp();
			int strideNode = idxr.strideColNode();

			// Danckwert boundary condition implementation
			//cache.boundary[0] = boundFunc(t, comp, para); // c_in -> inlet DOF idas suggestion
			disc.boundary[1] = c[disc.nCol * strideCell - strideNode]; // c_r outlet
			disc.boundary[2] = -disc.g[0]; // S_l inlet
			disc.boundary[3] = -disc.g[disc.nCol * strideCell - strideNode]; // S_r outlet
		}

		double inStream(double t, const Discretization& _disc, const unsigned int comp, const unsigned int secIdx) {

			double in = 0.0;
			for (int i = 0; i < 4; i++) {
				in = in * t + _disc.coefs[comp](secIdx, 3 - i);
			}

			return in;
		}

		void consistentInitialization(Discretization& disc, double t, double* const yPtr, double* const ypPtr, unsigned int secIdx) {

			unsigned int DOFs = disc.nPoints * disc.nComp;
			unsigned int nPoints = disc.nPoints;

			const double* cPtr = yPtr + disc.nComp;

			Eigen::Map<Eigen::VectorXd> yp(ypPtr, disc.nComp + (2 * DOFs));

			// Calculate solid phase RHS
			calcRHSq_DG(yPtr, ypPtr, disc);

			yp.segment(0, disc.nComp).setZero(); // TODO: time derivative of inlet function, 0.0 for inlet pulse for now

			for (int comp = 0; comp < disc.nComp; comp++) {
				// Convection dispersion RHS for one component
				ConvDisp_DG(cPtr, ypPtr + disc.nComp, disc, t, comp, secIdx);
				// Residual plugged into yp
				yp.segment(disc.nComp + comp * nPoints, nPoints) = -yp.segment(disc.nComp, nPoints)
															     - yp.segment(disc.nComp + disc.nComp * nPoints + comp * nPoints, nPoints) * ((1 - disc.porosity) / disc.porosity);

				// Residual plugged into yp
				if (!disc.isKinetic[comp]) {
					yp.segment(disc.nComp + DOFs + comp* disc.nPoints, disc.nPoints).setZero();
				}
				else {
					// residual already stored in yp at q
				}
			}
		}

		// ==========================================================================================================================================================  //
		// ==========================================================================================================================================================  //
		// ========================================						DG Jacobian							======================================================  //
		// ==========================================================================================================================================================  //
		// ==========================================================================================================================================================  //


		// @TODO: for sparse jacobian
		int calcNNZ() {
			if (_disc.modal)
				return 0;
			else
				return 0;
		}

		int calcStaticAnaJacobian(int secIdx) {

			// reset
			_jac.setZero();

			// DG convection dispersion Jacobian
			if (_disc.modal)
				calcStaticAnaModalJacobian(secIdx);
			else
				calcStaticAnaNodalJacobian(secIdx);

			// inlet DOFs Jacobian ( forward flow! )
			_jac.block(0, 0, _disc.nComp, _disc.nComp) = MatrixXd::Identity(_disc.nComp, _disc.nComp);//_jacInlet;

			// isotherm Jacobian
			calcIsothermJacobian();

			return 0;
		}

		int calcStaticAnaNodalJacobian(int secIdx) {

			unsigned int nNodes = _disc.nNodes;
			unsigned int nCells = _disc.nCol;
			unsigned int polyDeg = _disc.polyDeg;
			unsigned int nComp = _disc.nComp;

			// @TODO: special cases?
			if (nCells < 3)
				throw std::invalid_argument("Nodal Jacobian special case for nCells < 3 not implemented (yet?)");

			// inlet DOFs -> separate block: _jacInlet
			//_jac.block(0, 0, nComp, nComp) = MatrixXd::Identity(nComp, nComp);

			/*			Define inner cell Convection and Dispersion blocks			*/

			// Convection block [ d RHS_conv / d c ], depends only on first entry of previous cell
			MatrixXd convBlock = MatrixXd::Zero(nNodes, nNodes + 1);
			convBlock.block(0, 1, nNodes, nNodes) = _disc.polyDerM;
			convBlock(0, 0) = -_disc.invWeights[0];
			convBlock(0, 1) += _disc.invWeights[0];
			convBlock *= 2 * _disc.velocity[secIdx] / _disc.deltaZ;

			_jac.block(nComp, 0, nNodes, 1) = convBlock.block(0, 0, nNodes, 1); // inlet DOF
			_jac.block(nComp, nComp, nNodes, nNodes) = convBlock.block(0, 1, nNodes, nNodes);
			for (int cell = 1; cell < nCells; cell++) {
				_jac.block(nComp + cell * nNodes, nComp - 1 + cell * nNodes, nNodes, nNodes + 1) = convBlock;
			}

			// auxiliary Block for [ d g(c) / d c ], needed in Dispersion block
			MatrixXd GBlock = MatrixXd::Zero(nNodes, nNodes + 2);
			GBlock.block(0, 1, nNodes, nNodes) = _disc.polyDerM;
			GBlock(0, 0) -= 0.5 * _disc.invWeights[0];
			GBlock(0, 1) += 0.5 * _disc.invWeights[0];
			GBlock(nNodes - 1, nNodes) -= 0.5 * _disc.invWeights[polyDeg];
			GBlock(nNodes - 1, nNodes + 1) += 0.5 * _disc.invWeights[polyDeg];
			GBlock *= 2 * std::sqrt(_disc.dispersion[secIdx]) / _disc.deltaZ;

			// auxiliary for sum_{k=0}^{N} D_i,k G_k,j for all i, j
			MatrixXd sum = MatrixXd::Zero(nNodes, nNodes + 2);
			sum = _disc.polyDerM * GBlock;

			// Dispersion block [ d RHS_disp / d c ], depends on whole previous and subsequent cell
			MatrixXd dispBlock = MatrixXd::Zero(nNodes, 3 * nNodes); //
			// NOTE: N = polyDeg
			//cell indices : 0	 , ..., nNodes - 1;	nNodes, ..., 2 * nNodes - 1;	2 * nNodes, ..., 3 * nNodes - 1
			//			j  : -N-1, ..., -1		  ; 0     , ..., N			   ;	N + 1, ..., 2N + 1
			dispBlock.block(0, nNodes - 1, nNodes, nNodes + 2) = sum;
			dispBlock(0, nNodes - 1) += -_disc.invWeights[0] * (-0.5 * GBlock(0, 0) + 0.5 * GBlock(nNodes - 1, nNodes)); // G_N,N		i=0, j=-1
			dispBlock(0, nNodes) += -_disc.invWeights[0] * (-0.5 * GBlock(0, 1) + 0.5 * GBlock(nNodes - 1, nNodes + 1)); // G_N,N+1	i=0, j=0
			dispBlock.block(0, nNodes + 1, 1, nNodes) += -_disc.invWeights[0] * (-0.5 * GBlock.block(0, 2, 1, nNodes)); // G_i,j		i=0, j=1,...,N+1
			dispBlock.block(0, 0, 1, nNodes - 1) += -_disc.invWeights[0] * (0.5 * GBlock.block(nNodes - 1, 1, 1, nNodes - 1)); // G_N,j+N+1		i=0, j=-N-1,...,-2
			dispBlock.block(nNodes - 1, nNodes - 1, 1, nNodes) += _disc.invWeights[nNodes - 1] * (-0.5 * GBlock.block(nNodes - 1, 0, 1, nNodes)); // G_i,j+N+1		i=N, j=--1,...,N-1
			dispBlock(nNodes - 1, 2 * nNodes - 1) += _disc.invWeights[nNodes - 1] * (-0.5 * GBlock(nNodes - 1, nNodes) + 0.5 * GBlock(0, 0)); // G_i,j		i=N, j=N
			dispBlock(nNodes - 1, 2 * nNodes) += _disc.invWeights[nNodes - 1] * (-0.5 * GBlock(nNodes - 1, nNodes + 1) + 0.5 * GBlock(0, 1)); // G_i,j		i=N, j=N+1
			dispBlock.block(nNodes - 1, 2 * nNodes + 1, 1, nNodes - 1) += _disc.invWeights[nNodes - 1] * (0.5 * GBlock.block(0, 2, 1, nNodes - 1)); // G_0,j-N-1		i=N, j=N+2,...,2N+1
			dispBlock *= 2 * std::sqrt(_disc.dispersion[secIdx]) / _disc.deltaZ;

			// insert Blocks to Jacobian inner cells
			for (int cell = 1; cell < nCells - 1; cell++) {
				_jac.block(nComp + cell * nNodes, nComp + (cell - 1) * nNodes, nNodes, 3 * nNodes) -= dispBlock;
			}

			/*			Define boundary cell Convection and Dispersion blocks			*/
			/* left cell */
			// adjust auxiliary Block [ d g(c) / d c ] for left boundary cell
			MatrixXd GBlockBound = GBlock;
			GBlockBound(0, 1) -= 0.5 * _disc.invWeights[0] * 2 * std::sqrt(_disc.dispersion[secIdx]) / _disc.deltaZ;
			sum = _disc.polyDerM * GBlockBound; // correct corresponding column in D*G

			// estimate dispersion block ( j < 0 not needed)
			dispBlock.setZero();
			dispBlock.block(0, nNodes - 1, nNodes, nNodes + 2) = sum;
			dispBlock.block(0, nNodes - 1, 1, nNodes + 2) += -_disc.invWeights[0] * (-GBlockBound.block(0, 0, 1, nNodes + 2)); // G_N,N		i=0, j=-1,...,N+1
			dispBlock.block(nNodes - 1, nNodes - 1, 1, nNodes) += _disc.invWeights[nNodes - 1] * (-0.5 * GBlockBound.block(nNodes - 1, 0, 1, nNodes)); // G_i,j+N+1		i=N, j=--1,...,N-1
			dispBlock(nNodes - 1, 2 * nNodes - 1) += _disc.invWeights[nNodes - 1] * (-0.5 * GBlockBound(nNodes - 1, nNodes) + 0.5 * GBlockBound(0, 0)); // G_i,j		i=N, j=N
			dispBlock(nNodes - 1, 2 * nNodes) += _disc.invWeights[nNodes - 1] * (-0.5 * GBlockBound(nNodes - 1, nNodes + 1) + 0.5 * GBlock(0, 1)); // G_i,j		i=N, j=N+1
			dispBlock.block(nNodes - 1, 2 * nNodes + 1, 1, nNodes - 1) += _disc.invWeights[nNodes - 1] * (0.5 * GBlock.block(0, 2, 1, nNodes - 1)); // G_0,j-N-1		i=N, j=N+2,...,2N+1
			dispBlock *= 2 * std::sqrt(_disc.dispersion[secIdx]) / _disc.deltaZ;
			// copy *-1 to Jacobian
			_jac.block(nComp, nComp, nNodes, 2 * nNodes) -= dispBlock.block(0, nNodes, nNodes, 2 * nNodes);

			/* right cell */
		   // adjust auxiliary Block [ d g(c) / d c ] for left boundary cell
			GBlockBound(0, 1) += 0.5 * _disc.invWeights[0] * 2 * std::sqrt(_disc.dispersion[secIdx]) / _disc.deltaZ; 	// reverse change from left boundary
			GBlockBound(nNodes - 1, nNodes) += 0.5 * _disc.invWeights[polyDeg] * 2 * std::sqrt(_disc.dispersion[secIdx]) / _disc.deltaZ;
			sum = _disc.polyDerM * GBlockBound; // correct corresponding column in D*G

			// estimate dispersion block (only estimation differences to inner cell at N = 0 and j > N not needed)
			dispBlock.block(0, nNodes - 1, nNodes, nNodes + 2) = sum;
			dispBlock(0, nNodes - 1) += -_disc.invWeights[0] * (-0.5 * GBlockBound(0, 0) + 0.5 * GBlock(nNodes - 1, nNodes)); // G_N,N		i=0, j=-1
			dispBlock(0, nNodes) += -_disc.invWeights[0] * (-0.5 * GBlockBound(0, 1) + 0.5 * GBlock(nNodes - 1, nNodes + 1)); // G_N,N+1	i=0, j=0
			dispBlock.block(0, nNodes + 1, 1, nNodes) += -_disc.invWeights[0] * (-0.5 * GBlockBound.block(0, 2, 1, nNodes)); // G_i,j		i=0, j=1,...,N+1
			dispBlock.block(0, 0, 1, nNodes - 1) += -_disc.invWeights[0] * (0.5 * GBlock.block(nNodes - 1, 1, 1, nNodes - 1)); // G_N,j+N+1		i=0, j=-N-1,...,-2
			dispBlock.block(nNodes - 1, nNodes - 1, 1, nNodes + 2) += _disc.invWeights[nNodes - 1] * (-GBlockBound.block(nNodes - 1, 0, 1, nNodes + 2)); // G_i,j+N+1		i=N, j=--1,...,N+1
			dispBlock *= 2 * std::sqrt(_disc.dispersion[secIdx]) / _disc.deltaZ;
			// copy *-1 to Jacobian
			_jac.block(nComp + (nCells - 1) * nNodes, nComp + (nCells - 1 - 1) * nNodes, nNodes, 2 * nNodes) -= dispBlock.block(0, 0, nNodes, 2 * nNodes);


			/*	Copy DG Block to to all components 	*/

			int compBlock = _disc.nPoints;
			for (int comp = 1; comp < nComp; comp++) {
				_jac.block(nComp + comp * compBlock, comp, compBlock, 1) = _jac.block(nComp, 0, compBlock, 1); // inlet DOF
				_jac.block(nComp + comp * compBlock, nComp + comp * compBlock, compBlock, compBlock) = _jac.block(nComp, nComp, compBlock, compBlock);
			}

			return 0;
		}

		int calcStaticAnaModalJacobian(int secIdx) {

			unsigned int nNodes = _disc.nNodes;
			unsigned int nCells = _disc.nCol;
			unsigned int polyDeg = _disc.polyDeg;
			unsigned int nComp = _disc.nComp;

			// @TODO: special cases?
			if (nCells < 5)
				throw std::invalid_argument("Modal Jacobian special case for nCells < 5 not implemented (yet?)");

			//// inlet DOFs -> separate Block _jacInlet
			//_jac.block(0, 0, nComp, nComp) = MatrixXd::Identity(nComp, nComp);

			/*======================================================*/
			/*			Define Convection Jacobian Block			*/
			/*======================================================*/

			// Convection block [ d RHS_conv / d c ], additionally depends on first entry of previous cell

			MatrixXd convBlock = MatrixXd::Zero(nNodes, nNodes + 1);
			convBlock.block(0, 0, nNodes, 1) += _disc.invMM.block(0, 0, nNodes, 1);
			convBlock.block(0, 1, nNodes, nNodes) -= _disc.polyDerM;
			convBlock.block(0, 1, nNodes, 1) -= _disc.invMM.block(0, 0, nNodes, 1);
			convBlock *= 2 * _disc.velocity[secIdx] / _disc.deltaZ;

			// insert convection Blocks to Jacobian inner cells
			_jac.block(nComp, 0, nNodes, 1) = -convBlock.block(0, 0, nNodes, 1); // inlet DOF
			_jac.block(nComp, nComp, nNodes, nNodes) = -convBlock.block(0, 1, nNodes, nNodes);
			for (int cell = 1; cell < nCells; cell++) {
				_jac.block(nComp + cell * nNodes, nComp - 1 + cell * nNodes, nNodes, nNodes + 1) = -convBlock;
			}


			/*======================================================*/
			/*			Define Dispersion Jacobian Block			*/
			/*======================================================*/

			/* Inner cells */

			// Auxiliary Block [ d g(c) / d c ], additionally depends on boundary entries of neighbouring cells
			MatrixXd gBlock = MatrixXd::Zero(nNodes, nNodes + 2);
			gBlock.block(0, 1, nNodes, nNodes) = _disc.polyDerM;
			gBlock.block(0, 0, nNodes, 1) -= 0.5 * _disc.invMM.block(0, 0, nNodes, 1);
			gBlock.block(0, 1, nNodes, 1) += 0.5 * _disc.invMM.block(0, 0, nNodes, 1);
			gBlock.block(0, nNodes, nNodes, 1) -= 0.5 * _disc.invMM.block(0, nNodes - 1, nNodes, 1);
			gBlock.block(0, nNodes + 1, nNodes, 1) += 0.5 * _disc.invMM.block(0, nNodes - 1, nNodes, 1);
			gBlock *= 2 * std::sqrt(_disc.dispersion[secIdx]) / _disc.deltaZ;

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
			// compute d g^* / d c
			gStarDC.block(0, nNodes, 1, nNodes + 2) += gBlock.block(0, 0, 1, nNodes + 2);
			gStarDC.block(0, 0, 1, nNodes + 2) += gBlock.block(/*0*/nNodes - 1, 0, 1, nNodes + 2);
			gStarDC.block(nNodes - 1, nNodes, 1, nNodes + 2) += gBlock.block(nNodes - 1, 0, 1, nNodes + 2);
			gStarDC.block(nNodes - 1, 2 * nNodes, 1, nNodes + 2) += gBlock.block(0, 0, 1, nNodes + 2);
			gStarDC *= 0.5;

			//  indices  dispBlock :   0	 ,   1   , ..., nNodes;	nNodes+1, ..., 2 * nNodes;	2*nNodes+1, ..., 3 * nNodes; 3*nNodes+1
			//	derivative index j  : -(N+1)-1, -(N+1),...,	 -1	  ;   0     , ...,		N	 ;	  N + 1	  , ..., 2N + 2    ; 2(N+1) +1
			dispBlock.block(0, nNodes, nNodes, nNodes + 2) += _disc.polyDerM * gBlock - _disc.invMM * B * gBlock;
			dispBlock += _disc.invMM * B * gStarDC;
			dispBlock *= 2 * std::sqrt(_disc.dispersion[secIdx]) / _disc.deltaZ;

			for (unsigned int cell = 2; cell < nCells - 2; cell++) {
				_jac.block(nComp + cell * nNodes, nComp + cell * nNodes - nNodes - 1, nNodes, 3 * nNodes + 2) -= dispBlock;
			}

			/*		boundary cell neighbours		*/

				// left boundary cell neighbour
				// boundary auxiliary block [ d g(c) / d c ]
			MatrixXd GBlockBound_l = MatrixXd::Zero(nNodes, nNodes + 2);
			GBlockBound_l.block(0, 1, nNodes, nNodes) += _disc.polyDerM;
			GBlockBound_l.block(0, nNodes, nNodes, 1) -= 0.5 * _disc.invMM.block(0, nNodes - 1, nNodes, 1);
			GBlockBound_l.block(0, nNodes + 1, nNodes, 1) += 0.5 * _disc.invMM.block(0, nNodes - 1, nNodes, 1);
			GBlockBound_l *= 2 * std::sqrt(_disc.dispersion[secIdx]) / _disc.deltaZ;

			gStarDC.setZero();
			gStarDC.block(0, nNodes, 1, nNodes + 2) += gBlock.block(0, 0, 1, nNodes + 2);
			gStarDC.block(0, 0, 1, nNodes + 2) += GBlockBound_l.block(nNodes - 1, 0, 1, nNodes + 2);
			gStarDC.block(nNodes - 1, nNodes, 1, nNodes + 2) += gBlock.block(nNodes - 1, 0, 1, nNodes + 2);
			gStarDC.block(nNodes - 1, 2 * nNodes, 1, nNodes + 2) += gBlock.block(0, 0, 1, nNodes + 2);
			gStarDC *= 0.5;

			dispBlock.setZero();
			dispBlock.block(0, nNodes, nNodes, nNodes + 2) += _disc.polyDerM * gBlock - _disc.invMM * B * gBlock;
			dispBlock += _disc.invMM * B * gStarDC;
			dispBlock *= 2 * std::sqrt(_disc.dispersion[secIdx]) / _disc.deltaZ;

			_jac.block(nComp + nNodes, nComp, nNodes, 3 * nNodes + 1) -= dispBlock.block(0, 1, nNodes, 3 * nNodes + 1);

			// right boundary cell neighbour
			// boundary auxiliary block [ d g(c) / d c ]
			MatrixXd GBlockBound_r = MatrixXd::Zero(nNodes, nNodes + 2);
			GBlockBound_r.block(0, 1, nNodes, nNodes) += _disc.polyDerM;
			GBlockBound_r.block(0, 0, nNodes, 1) -= 0.5 * _disc.invMM.block(0, 0, nNodes, 1);
			GBlockBound_r.block(0, 1, nNodes, 1) += 0.5 * _disc.invMM.block(0, 0, nNodes, 1);
			GBlockBound_r *= 2 * std::sqrt(_disc.dispersion[secIdx]) / _disc.deltaZ;

			gStarDC.setZero();
			gStarDC.block(0, nNodes, 1, nNodes + 2) += gBlock.block(0, 0, 1, nNodes + 2);
			gStarDC.block(0, 0, 1, nNodes + 2) += gBlock.block(nNodes - 1, 0, 1, nNodes + 2);
			gStarDC.block(nNodes - 1, nNodes, 1, nNodes + 2) += gBlock.block(nNodes - 1, 0, 1, nNodes + 2);
			gStarDC.block(nNodes - 1, 2 * nNodes, 1, nNodes + 2) += GBlockBound_r.block(0, 0, 1, nNodes + 2);
			gStarDC *= 0.5;

			dispBlock.setZero();
			dispBlock.block(0, nNodes, nNodes, nNodes + 2) += _disc.polyDerM * gBlock - _disc.invMM * B * gBlock;
			dispBlock += _disc.invMM * B * gStarDC;
			dispBlock *= 2 * std::sqrt(_disc.dispersion[secIdx]) / _disc.deltaZ;

			_jac.block(nComp + (nCells - 2) * nNodes, nComp + (nCells - 2) * nNodes - nNodes - 1, nNodes, 3 * nNodes + 2) -= dispBlock;

			/*			boundary cells			*/

			// left boundary cell
			dispBlock.setZero();
			gStarDC.setZero();
			gStarDC.block(nNodes - 1, nNodes, 1, nNodes + 2) += GBlockBound_l.block(nNodes - 1, 0, 1, nNodes + 2);
			gStarDC.block(nNodes - 1, 2 * nNodes, 1, nNodes + 2) += gBlock.block(0, 0, 1, nNodes + 2);
			gStarDC *= 0.5;
			dispBlock.block(0, nNodes, nNodes, nNodes + 2) += _disc.polyDerM * GBlockBound_l - _disc.invMM * B * GBlockBound_l;
			dispBlock.block(0, nNodes + 1, nNodes, 2 * nNodes + 1) += _disc.invMM * B * gStarDC.block(0, nNodes + 1, nNodes, 2 * nNodes + 1);
			dispBlock *= 2 * std::sqrt(_disc.dispersion[secIdx]) / _disc.deltaZ;

			_jac.block(nComp, nComp, nNodes, 2 * nNodes + 1) -= dispBlock.block(0, nNodes + 1, nNodes, 2 * nNodes + 1);

			// right boundary cell
			dispBlock.setZero();
			gStarDC.setZero();
			gStarDC.block(0, nNodes, 1, nNodes + 2) += GBlockBound_r.block(0, 0, 1, nNodes + 2);
			gStarDC.block(0, 0, 1, nNodes + 2) += gBlock.block(nNodes - 1, 0, 1, nNodes + 2);
			gStarDC *= 0.5;
			dispBlock.block(0, nNodes, nNodes, nNodes + 2) += _disc.polyDerM * GBlockBound_r - _disc.invMM * B * GBlockBound_r;
			dispBlock += _disc.invMM * B * gStarDC;
			dispBlock *= 2 * std::sqrt(_disc.dispersion[secIdx]) / _disc.deltaZ;

			_jac.block(nComp + (nCells - 1) * nNodes, nComp + (nCells - 1) * nNodes - nNodes - 1, nNodes, 2 * nNodes + 1)
				-= dispBlock.block(0, 0, nNodes, 2 * nNodes + 1);

			/*======================================================*/
			/*			Copy DG Block To All Components				*/
			/*======================================================*/

			int compBlock = _disc.nPoints;
			for (int comp = 1; comp < nComp; comp++) {
				_jac.block(nComp + comp * compBlock, comp, compBlock, 1) = _jac.block(nComp, 0, compBlock, 1); // inlet DOF
				_jac.block(nComp + comp * compBlock, nComp + comp * compBlock, compBlock, compBlock) = _jac.block(nComp, nComp, compBlock, compBlock);
			}

		}

		/* TODO: add to binding model implementation ? // no access to binding model ka, kd */
		int calcIsothermJacobian() {

			int compBlock = _disc.nPoints;
			int nComp = _disc.nComp;
			int DGBlock = compBlock * nComp;

			for (int comp = 0; comp < nComp; comp++) {
				if (_disc.isotherm == "LINEAR") {
					
					_jac.block(nComp + DGBlock + comp * compBlock, nComp + comp * compBlock, compBlock, compBlock)
						= -MatrixXd::Identity(compBlock, compBlock) * _disc.adsorption[comp];
					_jac.block(nComp + DGBlock + comp * compBlock, nComp + DGBlock + comp * compBlock, compBlock, compBlock)
						= MatrixXd::Identity(compBlock, compBlock) * ((_disc.isKinetic[comp]) ? 1.0 : _disc.desorption[comp]);
				}
				else {
					throw std::invalid_argument("isotherm not implemented yet!");
				}
			}

			return 0;
		}


		/* State derivative Jacobian*/

		void calcStatederJacobian(double c_j) {

			int nNodes = _disc.nNodes;
			int nComp = _disc.nComp;
			int nCells = _disc.nCol;
			int polyDeg = _disc.polyDeg;
			int DOFs = nNodes * nComp * nCells;

			// =================================================================================================== //
			//	 State derivative Jacobian: d Residual / d y_dot												   //
			// =================================================================================================== //
			// state derivative Jacobian -> identity matrix blocks only !
			int DGBlock = nComp * nCells * nNodes;
			for (int comp = 0; comp < nComp; comp++) {
				int strideComp = comp * (nCells * nNodes);
				//cell block of d rhs/ d c_t
				_jacDisc.block(nComp + strideComp, nComp + strideComp, nNodes * nCells, nNodes * nCells)
					+= c_j * MatrixXd::Identity(nNodes * nCells, nNodes * nCells);

				// cell block of d rhs/ d q_t
				_jacDisc.block(nComp + strideComp, nComp + DGBlock + strideComp, nNodes * nCells, nNodes * nCells)
					+= c_j * ((1 - _disc.porosity) / _disc.porosity) * MatrixXd::Identity(nNodes * nCells, nNodes * nCells);

				// cell block of d isotherm / d q_t
				if (_disc.isKinetic[comp]) {
					_jacDisc.block(nComp + DGBlock + strideComp, nComp + DGBlock + strideComp, nNodes * nCells, nNodes * nCells)
						+= c_j * MatrixXd::Identity(nNodes * nCells, nNodes * nCells);
				}
			}
		}


};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_LUMPEDRATEMODELWITHOUTPORESDG_HPP_
