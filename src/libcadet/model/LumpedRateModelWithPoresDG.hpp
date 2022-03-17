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
 * Defines the lumped rate model with pores (LRMP).
 */

#ifndef LIBCADET_LUMPEDRATEMODELWITHPORESDG_HPP_
#define LIBCADET_LUMPEDRATEMODELWITHPORESDG_HPP_

#include "UnitOperationBase.hpp"
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

#include "C:\Users\jmbr\Cadet\libs\eigen-3.4.0\Eigen\Dense.hpp"	// use LA lib Eigen for Matrix operations
#include "C:\Users\jmbr\Cadet\libs\eigen-3.4.0\Eigen\Sparse.hpp"
#include <array>
#include <vector>
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <math.h>
#endif

#include "Benchmark.hpp"

namespace cadet
{

namespace model
{

class IDynamicReactionModel;

/**
 * @brief Lumped rate model of liquid column chromatography with pores
 * @details See @cite Guiochon2006, @cite Gu1995, @cite Felinger2004
 *
 * @f[\begin{align}
	\frac{\partial c_i}{\partial t} &= - u \frac{\partial c_i}{\partial z} + D_{\text{ax},i} \frac{\partial^2 c_i}{\partial z^2} - \frac{1 - \varepsilon_c}{\varepsilon_c} \frac{3 k_{f,i}}{r_p} j_{f,i} \\
	\frac{\partial c_{p,i}}{\partial t} + \frac{1 - \varepsilon_p}{\varepsilon_p} \frac{\partial q_{i}}{\partial t} &= \frac{3 k_{f,i}}{\varepsilon_p r_p} j_{f,i} \\
	a \frac{\partial q_i}{\partial t} &= f_{\text{iso}}(c_p, q)
\end{align} @f]
@f[ \begin{align}
	j_{f,i} = c_i - c_{p,i}
\end{align} @f]
 * Danckwerts boundary conditions (see @cite Danckwerts1953)
@f[ \begin{align}
u c_{\text{in},i}(t) &= u c_i(t,0) - D_{\text{ax},i} \frac{\partial c_i}{\partial z}(t,0) \\
\frac{\partial c_i}{\partial z}(t,L) &= 0
\end{align} @f]
 * Methods are described in @cite VonLieres2010a (WENO, linear solver), @cite Puttmann2013 @cite Puttmann2016 (forward sensitivities, AD, band compression)
 */
class LumpedRateModelWithPoresDG : public UnitOperationBase
{
public:

	LumpedRateModelWithPoresDG(UnitOpIdx unitOpIdx);
	virtual ~LumpedRateModelWithPoresDG() CADET_NOEXCEPT;

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

	static const char* identifier() { return "LUMPED_RATE_MODEL_WITH_PORES_DG"; }
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
			_timerJacobianPar.totalElapsedTime(),
			_timerConsistentInit.totalElapsedTime(),
			_timerConsistentInitPar.totalElapsedTime(),
			_timerLinearSolve.totalElapsedTime(),
			_timerFactorize.totalElapsedTime(),
			_timerFactorizePar.totalElapsedTime(),
			_timerMatVec.totalElapsedTime(),
			_timerGmres.totalElapsedTime()
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
			"Gmres"
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
	int residualBulk(double t, unsigned int secIdx, StateType const* yBase, double const* yDotBase, ResidualType* resBase, util::ThreadLocalStorage& threadLocalMem);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualParticle(double t, unsigned int parType, unsigned int colCell, unsigned int secIdx, StateType const* y, double const* yDot, ResidualType* res, util::ThreadLocalStorage& threadLocalMem);

	template <typename StateType, typename ResidualType, typename ParamType>
	int residualFlux(double t, unsigned int secIdx, StateType const* y, double const* yDot, ResidualType* res);

	void assembleOffdiagJac(double t, unsigned int secIdx);
	void extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset);

	int schurComplementMatrixVector(double const* x, double* z) const;
	void assembleDiscretizedJacobianParticleBlock(unsigned int type, double alpha, const Indexer& idxr);

	void addTimeDerivativeToJacobianParticleBlock(linalg::FactorizableBandMatrix::RowIterator& jac, const Indexer& idxr, double alpha, unsigned int parType);
	void solveForFluxes(double* const vecState, const Indexer& idxr);

	unsigned int numAdDirsForJacobian() const CADET_NOEXCEPT;

	int multiplexInitialConditions(const cadet::ParameterId& pId, unsigned int adDirection, double adValue);
	int multiplexInitialConditions(const cadet::ParameterId& pId, double val, bool checkSens);

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	void checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const;
#endif

	class Discretization // @TODO: separate convDisp DGoperator class
	{
	public:
		unsigned int nComp; //!< Number of components
		unsigned int nCol; //!< Number of column cells
		unsigned int polyDeg; //!< polynomial degree
		unsigned int nNodes; //!< Number of nodes per cell
		unsigned int nPoints; //!< Number of discrete Points
		unsigned int nParType; //!< Number of particle types
		unsigned int* parTypeOffset; //!< Array with offsets (in particle block) to particle type, additional last element contains total number of particle DOFs
		unsigned int* nBound; //!< Array with number of bound states for each component and particle type (particle type major ordering)
		unsigned int* boundOffset; //!< Array with offset to the first bound state of each component in the solid phase (particle type major ordering)
		unsigned int* strideBound; //!< Total number of bound states for each particle type, additional last element contains total number of bound states for all types
		unsigned int* nBoundBeforeType; //!< Array with number of bound states before a particle type (cumulative sum of strideBound)

		//////////////////////		DG specifics		////////////////////////////////////////////////////#
		//
		//// NOTE: no different Riemann solvers or boundary conditions

		double deltaZ;
		unsigned int modal;	//!< bool switch: 1 for modal basis, 0 for nodal basis
		Eigen::VectorXd nodes; //!< Array with positions of nodes in reference element
		Eigen::MatrixXd polyDerM; //!< Array with polynomial derivative Matrix
		Eigen::VectorXd invWeights; //!< Array with weights for numerical quadrature of size nNodes
		Eigen::MatrixXd invMM; //!< dense !INVERSE! mass matrix for modal (exact) integration

		// vgl. convDispOperator
		Eigen::VectorXd dispersion; //!< Column dispersion (may be section and component dependent)
		bool _dispersionCompIndep; //!< Determines whether dispersion is component independent
		double velocity; //!< Interstitial velocity (may be section dependent) \f$ u \f$
		double curVelocity;
		double crossSection; //!< Cross section area 
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
		void barycentricWeights(Eigen::VectorXd& baryWeights) {
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
			Eigen::VectorXd baryWeights = Eigen::VectorXd::Ones(polyDeg + 1u);
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
		void jacobiPolynomial(const double alpha, const double beta, const int N, Eigen::MatrixXd& P, int index) {
			// factor needed to normalize the Jacobi polynomial using the gamma function
			double gamma0 = std::pow(2.0, alpha + beta + 1.0) / (alpha + beta + 1.0) * std::tgamma(alpha + 1.0) * std::tgamma(beta + 1) / std::tgamma(alpha + beta + 1);
			Eigen::MatrixXd PL(N + 1, nodes.size());
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
		Eigen::MatrixXd getVandermonde_JACOBI() {

			Eigen::MatrixXd V(nodes.size(), nodes.size());

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
		void GradVandermonde(Eigen::MatrixXd& V_x) {
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
			Eigen::MatrixXd V_xT = polyDerM.transpose();

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

	parts::ConvectionDispersionOperator _convDispOp; //!< Convection dispersion operator for interstitial volume transport
	IDynamicReactionModel* _dynReactionBulk; //!< Dynamic reactions in the bulk volume

	std::vector<linalg::BandMatrix> _jacP; //!< Particle jacobian diagonal blocks (all of them for each particle type)
	std::vector<linalg::FactorizableBandMatrix> _jacPdisc; //!< Particle jacobian diagonal blocks (all of them for each particle type) with time derivatives from BDF method

	linalg::DoubleSparseMatrix _jacCF; //!< Jacobian block connecting interstitial states and fluxes (interstitial transport equation)
	linalg::DoubleSparseMatrix _jacFC; //!< Jacobian block connecting fluxes and interstitial states (flux equation)
	std::vector<linalg::DoubleSparseMatrix> _jacPF; //!< Jacobian blocks connecting particle states and fluxes (particle transport boundary condition)
	std::vector<linalg::DoubleSparseMatrix> _jacFP; //!< Jacobian blocks connecting fluxes and particle states (flux equation)

	linalg::DoubleSparseMatrix _jacInlet; //!< Jacobian inlet DOF block matrix connects inlet DOFs to first bulk cells

	active _colPorosity; //!< Column porosity (external porosity) \f$ \varepsilon_c \f$
	std::vector<double> _parGeomSurfToVol; //!< Particle surface to volume ratio factor (i.e., 3.0 for spherical, 2.0 for cylindrical, 1.0 for hexahedral)
	std::vector<active> _parRadius; //!< Particle radius \f$ r_p \f$
	bool _singleParRadius;
	std::vector<active> _parPorosity; //!< Particle porosity (internal porosity) \f$ \varepsilon_p \f$
	bool _singleParPorosity;
	std::vector<active> _parTypeVolFrac; //!< Volume fraction of each particle type

	// Vectorial parameters
	std::vector<active> _filmDiffusion; //!< Film diffusion coefficient \f$ k_f \f$
	MultiplexMode _filmDiffusionMode;
	std::vector<active> _poreAccessFactor; //!< Pore accessibility factor \f$ F_{\text{acc}} \f$
	MultiplexMode _poreAccessFactorMode;

	bool _axiallyConstantParTypeVolFrac; //!< Determines whether particle type volume fraction is homogeneous across axial coordinate
	bool _analyticJac; //!< Determines whether AD or analytic Jacobians are used
	unsigned int _jacobianAdDirs; //!< Number of AD seed vectors required for Jacobian computation

	bool _factorizeJacobian; //!< Determines whether the Jacobian needs to be factorized
	double* _tempState; //!< Temporary storage with the size of the state vector or larger if binding models require it
	linalg::Gmres _gmres; //!< GMRES algorithm for the Schur-complement in linearSolve()
	double _schurSafety; //!< Safety factor for Schur-complement solution

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
		friend int schurComplementMultiplierLRMPoresDG(void* userData, double const* x, double* z);

	class Indexer
	{
	public:
		Indexer(const Discretization& disc) : _disc(disc) { }

		// Strides
		inline int strideColCell() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp); }
		inline int strideColComp() const CADET_NOEXCEPT { return 1; }

		inline int strideParComp() const CADET_NOEXCEPT { return 1; }
		inline int strideParLiquid() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp); }
		inline int strideParBound(int parType) const CADET_NOEXCEPT { return static_cast<int>(_disc.strideBound[parType]); }
		inline int strideParBlock(int parType) const CADET_NOEXCEPT { return strideParLiquid() + strideParBound(parType); }

		inline int strideFluxCell() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp) * static_cast<int>(_disc.nParType); }
		inline int strideFluxParType() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp); }
		inline int strideFluxComp() const CADET_NOEXCEPT { return 1; }

		// Offsets
		inline int offsetC() const CADET_NOEXCEPT { return _disc.nComp; }
		inline int offsetCp() const CADET_NOEXCEPT { return _disc.nComp * _disc.nCol + offsetC(); }
		inline int offsetCp(ParticleTypeIndex pti) const CADET_NOEXCEPT { return offsetCp() + _disc.parTypeOffset[pti.value]; }
		inline int offsetCp(ParticleTypeIndex pti, ParticleIndex pi) const CADET_NOEXCEPT { return offsetCp() + _disc.parTypeOffset[pti.value] + strideParBlock(pti.value) * pi.value; }
		inline int offsetJf() const CADET_NOEXCEPT { return offsetCp() + _disc.parTypeOffset[_disc.nParType]; }
		inline int offsetJf(ParticleTypeIndex pti) const CADET_NOEXCEPT { return offsetJf() + pti.value * _disc.nCol * _disc.nComp; }
		inline int offsetBoundComp(ParticleTypeIndex pti, ComponentIndex comp) const CADET_NOEXCEPT { return _disc.boundOffset[pti.value * _disc.nComp + comp.value]; }

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

		Exporter(const Discretization& disc, const LumpedRateModelWithPoresDG& model, double const* data) : _disc(disc), _idx(disc), _model(model), _data(data) { }
		Exporter(const Discretization&& disc, const LumpedRateModelWithPoresDG& model, double const* data) = delete;

		virtual bool hasParticleFlux() const CADET_NOEXCEPT { return true; }
		virtual bool hasParticleMobilePhase() const CADET_NOEXCEPT { return true; }
		virtual bool hasSolidPhase() const CADET_NOEXCEPT { return _disc.strideBound[_disc.nParType] > 0; }
		virtual bool hasVolume() const CADET_NOEXCEPT { return false; }

		virtual unsigned int numComponents() const CADET_NOEXCEPT { return _disc.nComp; }
		virtual unsigned int numAxialCells() const CADET_NOEXCEPT { return _disc.nCol; }
		virtual unsigned int numRadialCells() const CADET_NOEXCEPT { return 1u; }
		virtual unsigned int numInletPorts() const CADET_NOEXCEPT { return 1; }
		virtual unsigned int numOutletPorts() const CADET_NOEXCEPT { return 1; }
		virtual unsigned int numParticleTypes() const CADET_NOEXCEPT { return _disc.nParType; }
		virtual unsigned int numParticleShells(unsigned int parType) const CADET_NOEXCEPT { return 1u; }
		virtual unsigned int numBoundStates(unsigned int parType) const CADET_NOEXCEPT { return _disc.strideBound[parType]; }
		virtual unsigned int numBulkDofs() const CADET_NOEXCEPT { return _disc.nComp * _disc.nCol; }
		virtual unsigned int numParticleMobilePhaseDofs(unsigned int parType) const CADET_NOEXCEPT { return _disc.nComp * _disc.nCol; }
		virtual unsigned int numSolidPhaseDofs(unsigned int parType) const CADET_NOEXCEPT { return _disc.strideBound[parType] * _disc.nCol; }
		virtual unsigned int numFluxDofs() const CADET_NOEXCEPT { return _disc.nComp * _disc.nCol * _disc.nParType; }
		virtual unsigned int numVolumeDofs() const CADET_NOEXCEPT { return 0; }

		virtual double const* concentration() const { return _idx.c(_data); }
		virtual double const* flux() const { return _idx.jf(_data); }
		virtual double const* particleMobilePhase(unsigned int parType) const { return _data + _idx.offsetCp(ParticleTypeIndex{ parType }); }
		virtual double const* solidPhase(unsigned int parType) const { return _data + _idx.offsetCp(ParticleTypeIndex{ parType }) + _idx.strideParLiquid(); }
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
		virtual unsigned int particleMobilePhaseStride(unsigned int parType) const { return _idx.strideParBlock(parType); }
		virtual unsigned int solidPhaseStride(unsigned int parType) const { return _idx.strideParBlock(parType); }

		virtual void axialCoordinates(double* coords) const
		{
			const double h = static_cast<double>(_model._convDispOp.columnLength()) / static_cast<double>(_disc.nCol);
			for (unsigned int i = 0; i < _disc.nCol; ++i)
				coords[i] = (i + 0.5) * h;
		}
		virtual void radialCoordinates(double* coords) const { }
		virtual void particleCoordinates(unsigned int parType, double* coords) const
		{
			coords[0] = static_cast<double>(_model._parRadius[parType]) * 0.5;
		}

	protected:
		const Discretization& _disc;
		const Indexer _idx;
		const LumpedRateModelWithPoresDG& _model;
		double const* const _data;

		const std::array<StateOrdering, 2> _concentrationOrdering = { { StateOrdering::AxialCell, StateOrdering::Component } };
		const std::array<StateOrdering, 3> _particleOrdering = { { StateOrdering::ParticleType, StateOrdering::AxialCell, StateOrdering::Component } };
		const std::array<StateOrdering, 4> _solidOrdering = { { StateOrdering::ParticleType, StateOrdering::AxialCell, StateOrdering::Component, StateOrdering::BoundState } };
		const std::array<StateOrdering, 3> _fluxOrdering = { { StateOrdering::ParticleType, StateOrdering::AxialCell, StateOrdering::Component } };
	};
};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_LUMPEDRATEMODELWITHPORESDG_HPP_