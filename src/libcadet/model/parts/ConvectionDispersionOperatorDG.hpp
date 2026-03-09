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
 * Defines the convection dispersion transport operator according to the discontinuous Galerkin discretization.
 */

#ifndef LIBCADET_CONVECTIONDISPERSIONOPERATORDG_HPP_
#define LIBCADET_CONVECTIONDISPERSIONOPERATORDG_HPP_

#include "ParamIdUtil.hpp"
#include "AutoDiff.hpp"
#include "Memory.hpp"
#include "SimulationTypes.hpp"
#include <ParamReaderHelper.hpp>
#include "linalg/BandedEigenSparseRowIterator.hpp"
#include "model/parts/DGToolbox.hpp"

#include <unordered_map>
#include <unordered_set>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

using namespace Eigen;

namespace cadet
{

	class IParameterProvider;
	class IConfigHelper;
	struct AdJacobianParams;
	struct SimulationTime;
	class IModel;

	namespace model
	{

		class IParameterParameterDependence;

		namespace parts
		{

			/**
			 * @brief Convection dispersion transport operator
			 * @details Implements the equation
			 *
			 * @f[\begin{align}
				\frac{\partial c_i}{\partial t} &= - u \frac{\partial c_i}{\partial z} + D_{\text{ax},i} \frac{\partial^2 c_i}{\partial z^2} \\
			\end{align} @f]
			 * with Danckwerts boundary conditions (see @cite Danckwerts1953)
			@f[ \begin{align}
			u c_{\text{in},i}(t) &= u c_i(t,0) - D_{\text{ax},i} \frac{\partial c_i}{\partial z}(t,0) \\
			\frac{\partial c_i}{\partial z}(t,L) &= 0
			\end{align} @f]
			 * Methods are described in @cite Breuer2023, and @cite Puttmann2013, @cite Puttmann2016 (forward sensitivities, AD, band compression)
			 *
			 * This class does not store the Jacobian. It only fills existing matrices given to its residual() functions.
			 * It assumes that there is no offset to the inlet in the local state vector and that the firsts element is placed
			 * directly after the inlet DOFs.
			 */
			class AxialConvectionDispersionOperatorBaseDG
			{
			public:

				AxialConvectionDispersionOperatorBaseDG();
				~AxialConvectionDispersionOperatorBaseDG() CADET_NOEXCEPT;

				void setFlowRates(const active& in, const active& out, const active& colPorosity) CADET_NOEXCEPT;

				bool configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, unsigned int nComp, int polynomial_integration_mode, unsigned int nelements, unsigned int polyDeg, unsigned int strideNode);
				bool configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters);
				bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, Eigen::MatrixXd& jacInlet);

				int residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, double* res, Eigen::SparseMatrix<double, Eigen::RowMajor>& jac);
				int residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, active* res, Eigen::SparseMatrix<double, Eigen::RowMajor>& jac);
				int residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, double* res, WithoutParamSensitivity);
				int residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, active* res, WithParamSensitivity);
				int residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, WithParamSensitivity);
				int residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, WithoutParamSensitivity);

				int calcTransportJacobian(Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian, Eigen::MatrixXd& jacInlet, const int bulkOffset = 0);
				typedef Eigen::Triplet<double> T;
				void convDispJacPattern(std::vector<T>& tripletList, const int bulkOffset = 0);
				unsigned int nConvDispEntries(bool pureNNZ = false);
				void multiplyWithDerivativeJacobian(const SimulationTime& simTime, double const* sDot, double* ret) const;
				void addTimeDerivativeToJacobian(double alpha, Eigen::SparseMatrix<double, Eigen::RowMajor>& jacDisc, unsigned int blockOffset = 0);

				inline const active& columnLength() const CADET_NOEXCEPT { return _colLength; }
				inline const active& currentVelocity(double pos) const CADET_NOEXCEPT { return _curVelocity; }
				inline bool forwardFlow() const CADET_NOEXCEPT { return _curVelocity >= 0.0; }

				inline double elemLeftBound(unsigned int idx) const CADET_NOEXCEPT { return idx * static_cast<double>(_deltaZ); }
				double relativeCoordinate(unsigned int idx) const CADET_NOEXCEPT
				{
					//const unsigned int element = floor(idx / _nNodes);
					//const unsigned int node = idx % _nNodes;
					// divide by column length to get relative position
					return (floor(idx / _nNodes) * static_cast<double>(_deltaZ) + 0.5 * static_cast<double>(_deltaZ) * (1.0 + _nodes[idx % _nNodes])) / static_cast<double>(_colLength);
				}

				inline const double* LGLnodes() const CADET_NOEXCEPT { return &_nodes[0]; }
				inline const active& currentVelocity() const CADET_NOEXCEPT { return _curVelocity; }
				inline const active* currentDispersion(const int secIdx) const CADET_NOEXCEPT { return getSectionDependentSlice(_colDispersion, _nComp, secIdx); }
				inline const bool dispersionCompIndep() const CADET_NOEXCEPT { return _dispersionCompIndep; }

				inline unsigned int nComp() const CADET_NOEXCEPT { return _nComp; }
				inline unsigned int nelements() const CADET_NOEXCEPT { return _nElem; }
				inline unsigned int nNodes() const CADET_NOEXCEPT { return _nNodes; }
				inline unsigned int nPoints() const CADET_NOEXCEPT { return _nPoints; }
				inline bool exactInt() const CADET_NOEXCEPT { return _exactInt; }

				// Indexer functionality:
				// Strides
				inline int strideColElement() const CADET_NOEXCEPT { return static_cast<int>(_strideElem); }
				inline int strideColNode() const CADET_NOEXCEPT { return static_cast<int>(_strideNode); }
				inline int strideColComp() const CADET_NOEXCEPT { return 1; }
				// Offsets
				inline int offsetC() const CADET_NOEXCEPT { return _nComp; }

				unsigned int jacobianLowerBandwidth() const CADET_NOEXCEPT;
				unsigned int jacobianUpperBandwidth() const CADET_NOEXCEPT;
				double inletJacobianFactor() const CADET_NOEXCEPT;
				unsigned int nJacEntries(bool pureNNZ);

				// @todo use more efficient seed vectors. currently, we treat the jacobian as banded, but the pattern is actually more sparse when multiple components are considered
				// (note that active type directions are limited)
				// We have different jacobian structure for exact integration and collocation DG scheme, i.e. we need different seed vectors
				// collocation DG: 2 * N_n * (strideNode) + 1 = total bandwidth (main diagonal entries maximally depend on the next and last N_n liquid phase entries of same component)
				//    ex. int. DG: 4 * N_n * (strideNode) + 1 = total bandwidth (main diagonal entries maximally depend on the next and last 2*N_n liquid phase entries of same component)
				int requiredADdirs() const CADET_NOEXCEPT { return (_exactInt) ? 4 * _nNodes * strideColNode() + 1 : 2 * _nNodes * strideColNode() + 1; }


				bool setParameter(const ParameterId& pId, double value);
				bool setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue);
				bool setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& id, double value);

			protected:

				template <typename StateType, typename ResidualType, typename ParamType, typename RowIteratorType, bool wantJac>
				int residualImpl(const IModel& model, double t, unsigned int secIdx, StateType const* y, double const* yDot, ResidualType* res, RowIteratorType jacBegin);

				// discretization parameters
				unsigned int _nComp; //!< Number of components
				bool _exactInt;	//!< specifies whether integrals are calculated exactly or approximated with LGL quadrature
				unsigned int _polyDeg; //!< DG discretization polynomial degree
				unsigned int _nElem; //!< Number of axial elements
				unsigned int _nNodes; //!< Number of nodes per element
				unsigned int _nPoints; //!< Number of axial discrete points

				unsigned int _strideNode; //!< Number of values between the same item in two adjacent nodes
				unsigned int _strideElem; //!< Number of values between the same item in two adjacent Elements

				// discretization toolbox and memory buffers
				active _deltaZ; //!< element spacing
				Eigen::VectorXd _nodes; //!< LGL nodes in [-1, 1]
				Eigen::MatrixXd _polyDerM; //!< Polynomial derivative Matrix
				Eigen::VectorXd _invWeights; //!< Inverse LGL quadrature weights -> diagonal (lumped) LGL mass matrix
				Eigen::MatrixXd _invMM; //!< Inverse (exact) mass matrix
				Eigen::MatrixXd* _DGjacAxDispBlocks; //!< Unique Jacobian blocks for axial dispersion
				Eigen::MatrixXd _DGjacAxConvBlock; //!< Unique Jacobian blocks for axial convection

				active* _auxState; //!< auxiliary variable
				active* _subsState; //!< auxiliary substitute
				Eigen::Vector<active, Eigen::Dynamic> _surfaceFlux; //!< stores the surface flux values
				Eigen::Vector<active, 4> _boundary; //!< stores the boundary values from Danckwert boundary conditions

				// Simulation parameters
				active _colLength; //!< Column length \f$ L \f$
				active _crossSection; //!< Cross section area 

				// Section dependent parameters
				std::vector<active> _colDispersion; //!< Column dispersion (may be section dependent) \f$ D_{\text{ax}} \f$
				std::vector<active> _velocity; //!< Interstitial velocity (may be section dependent) \f$ u \f$
				active _curVelocity; //!< Current interstitial velocity \f$ u \f$ in this time section
				int _dir; //!< Current flow direction in this time section

				// needed?
				int _curSection; //!< current section index
				bool _newStaticJac; //!< determines wether static analytical jacobian needs to be computed (every section)

				bool _dispersionCompIndep; //!< Determines whether dispersion is component independent
				IParameterParameterDependence* _dispersionDep;

				/* ===================================================================================
				 *  Functions to calculate Jacobian blocks
				 * =================================================================================== */

				/**
				 * @brief calculates the convection part of the DG jacobian
				 */
				Eigen::MatrixXd DGjacobianConvBlock()
				{
					// Convection block [ d RHS_conv / d c ], additionally depends on upwind flux part from corresponding neighbour element
					Eigen::MatrixXd convBlock = Eigen::MatrixXd::Zero(_nNodes, _nNodes + 1);

					if (_curVelocity >= 0.0) { // forward flow -> Convection block additionally depends on last entry of previous element
						convBlock.block(0, 1, _nNodes, _nNodes) -= _polyDerM;

						if (_exactInt) {
							convBlock.block(0, 0, _nNodes, 1) += _invMM.block(0, 0, _nNodes, 1);
							convBlock.block(0, 1, _nNodes, 1) -= _invMM.block(0, 0, _nNodes, 1);
						}
						else {
							convBlock(0, 0) += _invWeights[0];
							convBlock(0, 1) -= _invWeights[0];
						}
					}
					else { // backward flow -> Convection block additionally depends on first entry of subsequent element
						convBlock.block(0, 0, _nNodes, _nNodes) -= _polyDerM;

						if (_exactInt) {
							convBlock.block(0, _nNodes - 1, _nNodes, 1) += _invMM.block(0, _nNodes - 1, _nNodes, 1);
							convBlock.block(0, _nNodes, _nNodes, 1) -= _invMM.block(0, _nNodes - 1, _nNodes, 1);
						}
						else {
							convBlock(_nNodes - 1, _nNodes - 1) += _invWeights[_nNodes - 1];
							convBlock(_nNodes - 1, _nNodes) -= _invWeights[_nNodes - 1];
						}
					}
					convBlock *= 2 / static_cast<double>(_deltaZ);

					return -convBlock; // *-1 for residual
				}

				/**
				 * @brief calculates the DG Jacobian auxiliary block
				 * @param [in] exInt true if exact integration DG scheme
				 * @param [in] elementIdx element index
				 */
				Eigen::MatrixXd getGBlock(unsigned int elementIdx)
				{
					// Auxiliary Block [ d g(c) / d c ], additionally depends on boundary entries of neighbouring elements
					Eigen::MatrixXd gBlock = Eigen::MatrixXd::Zero(_nNodes, _nNodes + 2);
					gBlock.block(0, 1, _nNodes, _nNodes) = _polyDerM;
					if (_exactInt) {
						if (elementIdx != 1 && elementIdx != _nElem) {
							gBlock.block(0, 0, _nNodes, 1) -= 0.5 * _invMM.block(0, 0, _nNodes, 1);
							gBlock.block(0, 1, _nNodes, 1) += 0.5 * _invMM.block(0, 0, _nNodes, 1);
							gBlock.block(0, _nNodes, _nNodes, 1) -= 0.5 * _invMM.block(0, _nNodes - 1, _nNodes, 1);
							gBlock.block(0, _nNodes + 1, _nNodes, 1) += 0.5 * _invMM.block(0, _nNodes - 1, _nNodes, 1);
						}
						else if (elementIdx == 1) { // left
							if (elementIdx == _nElem)
								return gBlock * 2 / static_cast<double>(_deltaZ);
							;
							gBlock.block(0, _nNodes, _nNodes, 1) -= 0.5 * _invMM.block(0, _nNodes - 1, _nNodes, 1);
							gBlock.block(0, _nNodes + 1, _nNodes, 1) += 0.5 * _invMM.block(0, _nNodes - 1, _nNodes, 1);
						}
						else if (elementIdx == _nElem) { // right
							gBlock.block(0, 0, _nNodes, 1) -= 0.5 * _invMM.block(0, 0, _nNodes, 1);
							gBlock.block(0, 1, _nNodes, 1) += 0.5 * _invMM.block(0, 0, _nNodes, 1);
						}
						else if (elementIdx == 0 || elementIdx == _nElem + 1) {
							gBlock.setZero();
						}
						gBlock *= 2 / static_cast<double>(_deltaZ);
					}
					else {
						if (elementIdx == 0 || elementIdx == _nElem + 1)
							return Eigen::MatrixXd::Zero(_nNodes, _nNodes + 2);

						gBlock(0, 0) -= 0.5 * _invWeights[0];
						gBlock(0, 1) += 0.5 * _invWeights[0];
						gBlock(_nNodes - 1, _nNodes) -= 0.5 * _invWeights[_nNodes - 1];
						gBlock(_nNodes - 1, _nNodes + 1) += 0.5 * _invWeights[_nNodes - 1];
						gBlock *= 2 / static_cast<double>(_deltaZ);

						if (elementIdx == 1) {
							// adjust auxiliary Block [ d g(c) / d c ] for left boundary element
							gBlock(0, 1) -= 0.5 * _invWeights[0] * 2 / static_cast<double>(_deltaZ);
							if (elementIdx == _nElem) { // adjust for special case one element
								gBlock(0, 0) += 0.5 * _invWeights[0] * 2 / static_cast<double>(_deltaZ);
								gBlock(_nNodes - 1, _nNodes + 1) -= 0.5 * _invWeights[_nNodes - 1] * 2 / static_cast<double>(_deltaZ);
								gBlock(_nNodes - 1, _nNodes) += 0.5 * _invWeights[_polyDeg] * 2 / static_cast<double>(_deltaZ);
							}
						}
						else if (elementIdx == _nElem) {
							// adjust auxiliary Block [ d g(c) / d c ] for right boundary element
							gBlock(_nNodes - 1, _nNodes) += 0.5 * _invWeights[_polyDeg] * 2 / static_cast<double>(_deltaZ);
						}
					}

					return gBlock;
				}
				/**
				 * @brief calculates the num. flux part of a dispersion DG Jacobian block
				 * @param [in] elementIdx element index
				 * @param [in] leftG left neighbour auxiliary block
				 * @param [in] middleG neighbour auxiliary block
				 * @param [in] rightG neighbour auxiliary block
				 */
				Eigen::MatrixXd auxBlockGstar(unsigned int elementIdx, Eigen::MatrixXd leftG, Eigen::MatrixXd middleG, Eigen::MatrixXd rightG)
				{
					// auxiliary block [ d g^* / d c ], depends on whole previous and subsequent element plus first entries of subsubsequent elements
					Eigen::MatrixXd gStarDC = Eigen::MatrixXd::Zero(_nNodes, 3 * _nNodes + 2);
					// NOTE: N = polyDeg
					// indices  gStarDC    :     0   ,   1   , ..., _nNodes; _nNodes+1, ..., 2 * _nNodes;	2*_nNodes+1, ..., 3 * _nNodes; 3*_nNodes+1
					// derivative index j  : -(N+1)-1, -(N+1),... ,  -1   ;   0     , ...,		N	 ;	  N + 1	  , ..., 2N + 2    ; 2(N+1) +1
					// auxiliary block [d g^* / d c]
					if (elementIdx != 1) {
						gStarDC.block(0, _nNodes, 1, _nNodes + 2) += middleG.block(0, 0, 1, _nNodes + 2);
						gStarDC.block(0, 0, 1, _nNodes + 2) += leftG.block(_nNodes - 1, 0, 1, _nNodes + 2);
					}
					if (elementIdx != _nElem) {
						gStarDC.block(_nNodes - 1, _nNodes, 1, _nNodes + 2) += middleG.block(_nNodes - 1, 0, 1, _nNodes + 2);
						gStarDC.block(_nNodes - 1, 2 * _nNodes, 1, _nNodes + 2) += rightG.block(0, 0, 1, _nNodes + 2);
					}
					gStarDC *= 0.5;

					return gStarDC;
				}

				Eigen::MatrixXd getBMatrix() {

					Eigen::MatrixXd B = Eigen::MatrixXd::Zero(_nNodes, _nNodes);
					B(0, 0) = -1.0;
					B(_nNodes - 1, _nNodes - 1) = 1.0;

					return B;
				}

				/**
				 * @brief calculates the dispersion part of the DG jacobian
				 * @param [in] exInt true if exact integration DG scheme
				 * @param [in] elementIdx element index
				 */
				Eigen::MatrixXd DGjacobianDispBlock(unsigned int elementIdx)
				{
					int offC = 0; // inlet DOFs not included in Jacobian

					Eigen::MatrixXd dispBlock;

					if (_exactInt) {

						// Inner dispersion block [ d RHS_disp / d c ], depends on whole previous and subsequent element plus first entries of subsubsequent elements
						dispBlock = Eigen::MatrixXd::Zero(_nNodes, 3 * _nNodes + 2);

						Eigen::MatrixXd B = getBMatrix(); // "Lifting" matrix
						Eigen::MatrixXd gBlock = getGBlock(elementIdx); // current element auxiliary block matrix
						Eigen::MatrixXd gStarDC = auxBlockGstar(elementIdx, getGBlock(elementIdx - 1), gBlock, getGBlock(elementIdx + 1)); // Numerical flux block

						//  indices  dispBlock :   0	 ,   1   , ..., _nNodes;	_nNodes+1, ..., 2 * _nNodes;	2*_nNodes+1, ..., 3 * _nNodes; 3*_nNodes+1
						//	derivative index j  : -(N+1)-1, -(N+1),...,	 -1	  ;   0     , ...,		N	 ;	  N + 1	  , ..., 2N + 2    ; 2(N+1) +1
						dispBlock.block(0, _nNodes, _nNodes, _nNodes + 2) += _polyDerM * gBlock - _invMM * B * gBlock;
						dispBlock += _invMM * B * gStarDC;
						dispBlock *= 2 / static_cast<double>(_deltaZ);
					}
					else { // inexact integration collocation DGSEM

						dispBlock = Eigen::MatrixXd::Zero(_nNodes, 3 * _nNodes);
						Eigen::MatrixXd GBlockLeft = getGBlock(elementIdx - 1);
						Eigen::MatrixXd GBlock = getGBlock(elementIdx);
						Eigen::MatrixXd GBlockRight = getGBlock(elementIdx + 1);

						// Dispersion block [ d RHS_disp / d c ], depends on whole previous and subsequent element
						// NOTE: N = polyDeg
						// element indices :  0  , ..., _nNodes - 1;	_nNodes, ..., 2 * _nNodes - 1;	2 * _nNodes, ..., 3 * _nNodes - 1
						//			 j  : -N-1, ..., -1		  ; 0     , ..., N			   ;	N + 1, ..., 2N + 1
						dispBlock.block(0, _nNodes - 1, _nNodes, _nNodes + 2) = _polyDerM * GBlock;

						if (elementIdx > 1) {
							dispBlock(0, _nNodes - 1) += -_invWeights[0] * (-0.5 * GBlock(0, 0) + 0.5 * GBlockLeft(_nNodes - 1, _nNodes)); // G_N,N		i=0, j=-1
							dispBlock(0, _nNodes) += -_invWeights[0] * (-0.5 * GBlock(0, 1) + 0.5 * GBlockLeft(_nNodes - 1, _nNodes + 1)); // G_N,N+1	i=0, j=0
							dispBlock.block(0, _nNodes + 1, 1, _nNodes) += -_invWeights[0] * (-0.5 * GBlock.block(0, 2, 1, _nNodes)); // G_i,j		i=0, j=1,...,N+1
							dispBlock.block(0, 0, 1, _nNodes - 1) += -_invWeights[0] * (0.5 * GBlockLeft.block(_nNodes - 1, 1, 1, _nNodes - 1)); // G_N,j+N+1		i=0, j=-N-1,...,-2
						}
						else if (elementIdx == 1) { // left boundary element
							dispBlock.block(0, _nNodes - 1, 1, _nNodes + 2) += -_invWeights[0] * (-GBlock.block(0, 0, 1, _nNodes + 2)); // G_N,N		i=0, j=-1,...,N+1
						}
						if (elementIdx < _nElem) {
							dispBlock.block(_nNodes - 1, _nNodes - 1, 1, _nNodes) += _invWeights[_nNodes - 1] * (-0.5 * GBlock.block(_nNodes - 1, 0, 1, _nNodes)); // G_i,j+N+1		i=N, j=-1,...,N-1
							dispBlock(_nNodes - 1, 2 * _nNodes - 1) += _invWeights[_nNodes - 1] * (-0.5 * GBlock(_nNodes - 1, _nNodes) + 0.5 * GBlockRight(0, 0)); // G_i,j		i=N, j=N
							dispBlock(_nNodes - 1, 2 * _nNodes) += _invWeights[_nNodes - 1] * (-0.5 * GBlock(_nNodes - 1, _nNodes + 1) + 0.5 * GBlockRight(0, 1)); // G_i,j		i=N, j=N+1
							dispBlock.block(_nNodes - 1, 2 * _nNodes + 1, 1, _nNodes - 1) += _invWeights[_nNodes - 1] * (0.5 * GBlockRight.block(0, 2, 1, _nNodes - 1)); // G_0,j-N-1		i=N, j=N+2,...,2N+1
						}
						else if (elementIdx == _nElem) { // right boundary element
							dispBlock.block(_nNodes - 1, _nNodes - 1, 1, _nNodes + 2) += _invWeights[_nNodes - 1] * (-GBlock.block(_nNodes - 1, 0, 1, _nNodes + 2)); // G_i,j+N+1		i=N, j=--1,...,N+1
						}

						dispBlock *= 2 / static_cast<double>(_deltaZ);
					}

					return -dispBlock; // *-1 for residual
				}

				/* ===================================================================================
				*   Residual functions
				* =================================================================================== */
				/**
				* @brief calculates the volume Integral of the auxiliary equation
				* @param [in] state state vector
				* @param [in, out] stateDer state derivative vector
				* @detail performs matrix-vector multiplication optimized depending on state types and stores the result in stateDer.
				*/
				template<typename StateType, typename ResidualType>
				void volumeIntegral(Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state, Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>>& stateDer) {

					for (unsigned int element = 0; element < _nElem; element++) {

						// exploit Eigen3 performance if no mixed scalar types
						if constexpr (std::is_same_v<StateType, double>) {
							if constexpr (std::is_same_v<ResidualType, double>) {
								stateDer.segment(element * _nNodes, _nNodes) -= _polyDerM * state.segment(element * _nNodes, _nNodes);
							}
							else {
								stateDer.segment(element * _nNodes, _nNodes) -= (_polyDerM * state.segment(element * _nNodes, _nNodes)).template cast<ResidualType>();
							}
						}
						else { // both active types
							// todo use custom (mixed scalar-type) matrix vector multiplication?
							stateDer.segment(element * _nNodes, _nNodes) -= _polyDerM.template cast<active>() * state.segment(element * _nNodes, _nNodes);
						}
					}
				}
				template<typename StateType, typename ResidualType>
				void volumeIntegral(Eigen::Map<Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state, Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>>& stateDer) {
					Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> state_const(state.data(), state.size(), InnerStride<Dynamic>(1));
					volumeIntegral(state_const, stateDer);
				}
				/*
				 * @brief calculates the interface fluxes h* of Convection Dispersion equation
				 */
				template<typename StateType, typename ParamType>
				void InterfaceFlux(const StateType* C, ParamType _dispersion) {

					StateType* g = reinterpret_cast<StateType*>(_auxState);

					// component-wise strides
					unsigned int strideNode = _strideNode;
					unsigned int strideelement = _nNodes * strideNode;
					unsigned int strideNode_g = 1u;
					unsigned int strideelement_g = _nNodes * strideNode_g;

					// Conv.Disp. flux: h* = h*_conv + h*_disp = v c_up + 0.5 sqrt(D_ax) (S_l + S_r)

					if (_curVelocity >= 0.0) { // forward flow (upwind num. flux)
						// calculate inner interface fluxes
						for (unsigned int element = 1; element < _nElem; element++) {
							// h* = h*_conv + h*_disp
							_surfaceFlux[element] // inner interfaces
								= _curVelocity * (C[element * strideelement - strideNode]) // left element (i.e. forward flow upwind)
								- 0.5 * (-2.0 / static_cast<ParamType>(_deltaZ)) * _dispersion *
								(g[element * strideelement_g - strideNode_g] // left element
									+ g[element * strideelement_g]); // right element
						}

						// boundary fluxes
						// inlet (left) boundary interface
						_surfaceFlux[0]
							= _curVelocity * _boundary[0];

						// outlet (right) boundary interface
						_surfaceFlux[_nElem]
							= _curVelocity * (C[_nElem * strideelement - strideNode])
							- 0.5 * (-2.0 / static_cast<ParamType>(_deltaZ)) * _dispersion *
							(g[_nElem * strideelement_g - strideNode_g] // last element last node
								+ _boundary[3]); // right boundary value g
					}
					else { // backward flow (upwind num. flux)
						// calculate inner interface fluxes
						for (unsigned int element = 1; element < _nElem; element++) {
							// h* = h*_conv + h*_disp
							_surfaceFlux[element] // inner interfaces
								= _curVelocity * (C[element * strideelement]) // right element (i.e. backward flow upwind)
								- 0.5 * (-2.0 / static_cast<ParamType>(_deltaZ)) * _dispersion *
								(g[element * strideelement_g - strideNode_g] // left element
									+ g[element * strideelement_g]); // right element
						}

						// boundary fluxes
						// inlet boundary interface
						_surfaceFlux[_nElem]
							= _curVelocity * _boundary[0];

						// outlet boundary interface
						_surfaceFlux[0]
							= _curVelocity * (C[0])
							- 0.5 * (-2.0 / static_cast<ParamType>(_deltaZ)) * _dispersion *
							(g[0] // first element first node
								+ _boundary[2]); // left boundary value g
					}
					// apply inverse mapping jacobian (reference space)
					_surfaceFlux *= -2.0 / static_cast<ParamType>(_deltaZ);
				}
				/**
				 * @brief calculates and fills the surface flux values for auxiliary equation
				 * @param [in] C bulk liquid phase
				 * @param [in] strideNode node stride w.r.t. C
				 * @param [in] strideelement element stride w.r.t. C
				 */
				template<typename StateType>
				void InterfaceFluxAuxiliary(const StateType* C, unsigned int strideNode, unsigned int strideelement) {

					// Auxiliary flux: c* = 0.5 (c_l + c_r)

					// calculate inner interface fluxes
					for (unsigned int element = 1; element < _nElem; element++) {
						_surfaceFlux[element] // left interfaces
							= 0.5 * (C[element * strideelement - strideNode] + // left node
								C[element * strideelement]); // right node
					}
					// calculate boundary interface fluxes

					_surfaceFlux[0] // left boundary interface
						= 0.5 * (C[0] + // boundary value
							C[0]); // first element first node

					_surfaceFlux[(_nElem)] // right boundary interface
						= 0.5 * (C[_nElem * strideelement - strideNode] + // last element last node
							C[_nElem * strideelement - strideNode]);// // boundary value
				}
				/**
				 * @brief calculates the string form surface Integral
				 * @param [in] state relevant state vector
				 * @param [in] stateDer state derivative vector the solution is added to
				 * @param [in] strideNode_state node stride w.r.t. state
				 * @param [in] strideelement_state element stride w.r.t. state
				 * @param [in] strideNode_stateDer node stride w.r.t. stateDer
				 * @param [in] strideelement_stateDer element stride w.r.t. stateDer
				 * @detail calculates stateDer = M^-1 * B * (state - state^*) and exploits LGL sparsity if applied
				 */
				template<typename StateType, typename ResidualType>
				void surfaceIntegral(const StateType* state, ResidualType* stateDer,
					unsigned int strideNode_state, unsigned int strideelement_state,
					unsigned int strideNode_stateDer, unsigned int strideelement_stateDer) {

					if (_exactInt) { // non-collocated integration -> dense mass matrix
						for (unsigned int element = 0; element < _nElem; element++) {
							// strong surface integral -> M^-1 B [state - state*]
							for (unsigned int Node = 0; Node < _nNodes; Node++) {
								stateDer[element * strideelement_stateDer + Node * strideNode_stateDer]
									-= static_cast<ResidualType>(_invMM(Node, 0) * (state[element * strideelement_state] - _surfaceFlux[element])
										- _invMM(Node, _polyDeg) * (state[element * strideelement_state + _polyDeg * strideNode_state] - _surfaceFlux[(element + 1)]));
							}
						}
					}
					else { // collocated numerical integration -> diagonal mass matrix
						for (unsigned int element = 0; element < _nElem; element++) {
							// strong surface integral -> M^-1 B [state - state*]
							stateDer[element * strideelement_stateDer] // first element, node
								-= static_cast<ResidualType>(_invWeights[0]
									* (state[element * strideelement_state] - _surfaceFlux(element)));

							stateDer[element * strideelement_stateDer + _polyDeg * strideNode_stateDer] // last element, node
								+= static_cast<ResidualType>(_invWeights[_polyDeg]
									* (state[element * strideelement_state + _polyDeg * strideNode_state] - _surfaceFlux(element + 1)));
						}
					}
				}
				/**
				 * @brief computes ghost nodes to implement boundary conditions
				 * @detail to implement Danckwert boundary conditions, we only need to set the solid wall BC values for auxiliary variable
				 */
				template<typename StateType>
				void calcBoundaryValues() {
					//cache.boundary[0] = c_in -> inlet DOF already set
					//_boundary[1] = (_velocity >= 0.0) ? C[_nPoints - 1] : C[0]; // c_r outlet not required in Danckwerts BC
					_boundary[2] = -reinterpret_cast<StateType*>(_auxState)[0]; // g_l left boundary (inlet/outlet for forward/backward flow)
					_boundary[3] = -reinterpret_cast<StateType*>(_auxState)[_nPoints - 1]; // g_r right boundary (outlet/inlet for forward/backward flow)
				}

				// ==========================================================================================================================================================  //
				// ========================================						DG Jacobian							=========================================================  //
				// ==========================================================================================================================================================  //

				/**
				 * @brief sets the sparsity pattern of the convection dispersion Jacobian for the collocation DG scheme
				 */
				int ConvDispCollocationPattern(std::vector<T>& tripletList, const int offC = 0) {

					/*======================================================*/
					/*			Define Convection Jacobian Block			*/
					/*======================================================*/

					// Convection block [ d RHS_conv / d c ], also depends on upwind entry

					if (_curVelocity >= 0.0) { // forward flow upwind entry -> last node of previous element
					// special inlet DOF treatment for inlet boundary element (first element)
						for (unsigned int comp = 0; comp < _nComp; comp++) {
							for (unsigned int i = 0; i < _nNodes; i++) {
								//tripletList.push_back(T(offC + comp * sComp + i * sNode, comp * sComp, 0.0)); // inlet DOFs not included in Jacobian
								for (unsigned int j = 1; j < _nNodes + 1; j++) {
									tripletList.push_back(T(offC + comp * strideColComp() + i * strideColNode(),
										offC + comp * strideColComp() + (j - 1) * strideColNode(),
										0.0));
								}
							}
						}
						for (unsigned int element = 1; element < _nElem; element++) {
							for (unsigned int comp = 0; comp < _nComp; comp++) {
								for (unsigned int i = 0; i < _nNodes; i++) {
									for (unsigned int j = 0; j < _nNodes + 1; j++) {
										// row: jump over inlet DOFs and previous elements, add component offset and go node strides from there for each convection block entry
										// col: jump over inlet DOFs and previous elements, go back one node, add component offset and go node strides from there for each convection block entry
										tripletList.push_back(T(offC + element * strideColElement() + comp * strideColComp() + i * strideColNode(),
											offC + element * strideColElement() - strideColNode() + comp * strideColComp() + j * strideColNode(),
											0.0));
									}
								}
							}
						}
					}
					else { // backward flow upwind entry -> first node of subsequent element
						// special inlet DOF treatment for inlet boundary element (last element)
						for (unsigned int comp = 0; comp < _nComp; comp++) {
							for (unsigned int i = 0; i < _nNodes; i++) {
								// inlet DOFs not included in Jacobian
								for (unsigned int j = 0; j < _nNodes; j++) {
									tripletList.push_back(T(offC + (_nElem - 1) * strideColElement() + comp * strideColComp() + i * strideColNode(),
										offC + (_nElem - 1) * strideColElement() + comp * strideColComp() + j * strideColNode(),
										0.0));
								}
							}
						}
						for (unsigned int element = 0; element < _nElem - 1u; element++) {
							for (unsigned int comp = 0; comp < _nComp; comp++) {
								for (unsigned int i = 0; i < _nNodes; i++) {
									for (unsigned int j = 0; j < _nNodes + 1; j++) {
										// row: jump over inlet DOFs and previous elements, add component offset and go node strides from there for each convection block entry
										// col: jump over inlet DOFs and previous elements, add component offset and go node strides from there for each convection block entry
										tripletList.push_back(T(offC + element * strideColElement() + comp * strideColComp() + i * strideColNode(),
											offC + element * strideColElement() + comp * strideColComp() + j * strideColNode(),
											0.0));
									}
								}
							}
						}
					}

					/*======================================================*/
					/*			Define Dispersion Jacobian Block			*/
					/*======================================================*/

					/*		Inner element dispersion blocks		*/

					// Dispersion block [ d RHS_disp / d c ], depends on whole previous and subsequent element

					// insert Blocks to Jacobian inner elements (only for _nElem >= 3)
					if (_nElem >= 3u) {
						for (unsigned int element = 1; element < _nElem - 1; element++) {
							for (unsigned int comp = 0; comp < _nComp; comp++) {
								for (unsigned int i = 0; i < _nNodes; i++) {
									for (unsigned int j = 0; j < 3 * _nNodes; j++) {
										// pattern is more sparse than a nNodes x 3*nNodes block.
										if ((j >= _nNodes - 1 && j <= 2 * _nNodes) ||
											(i == 0 && j <= 2 * _nNodes) ||
											(i == _nNodes - 1 && j >= _nNodes - 1))
											// row: jump over inlet DOFs and previous elements, add component offset and go node strides from there for each entry
											// col: jump over inlet DOFs and previous elements, go back one element, add component offset and go node strides from there for each entry
											tripletList.push_back(T(offC + element * strideColElement() + comp * strideColComp() + i * strideColNode(),
												offC + (element - 1) * strideColElement() + comp * strideColComp() + j * strideColNode(),
												0.0));
									}
								}
							}
						}
					}

					/*				Boundary element Dispersion blocks			*/

					if (_nElem != 1) { // Note: special case _nElem = 1 already set by advection block
						for (unsigned int comp = 0; comp < _nComp; comp++) {
							for (unsigned int i = 0; i < _nNodes; i++) {
								for (unsigned int j = _nNodes; j < 3 * _nNodes; j++) {
									// pattern is more sparse than a nNodes x 2*nNodes block.
									if ((j >= _nNodes - 1 && j <= 2 * _nNodes) ||
										(i == 0 && j <= 2 * _nNodes) ||
										(i == _nNodes - 1 && j >= _nNodes - 1))
										tripletList.push_back(T(offC + comp * strideColComp() + i * strideColNode(),
											offC + comp * strideColComp() + (j - _nNodes) * strideColNode(),
											0.0));
								}
							}
						}

						for (unsigned int comp = 0; comp < _nComp; comp++) {
							for (unsigned int i = 0; i < _nNodes; i++) {
								for (unsigned int j = 0; j < 2 * _nNodes; j++) {
									// pattern is more sparse than a nNodes x 2*nNodes block.
									if ((j >= _nNodes - 1 && j <= 2 * _nNodes) ||
										(i == 0 && j <= 2 * _nNodes) ||
										(i == _nNodes - 1 && j >= _nNodes - 1))
										tripletList.push_back(T(offC + (_nElem - 1) * strideColElement() + comp * strideColComp() + i * strideColNode(),
											offC + (_nElem - 1 - 1) * strideColElement() + comp * strideColComp() + j * strideColNode(),
											0.0));
								}
							}
						}
					}

					return 0;
				}

				/**
				* @brief sets the sparsity pattern of the convection dispersion Jacobian for the exact integration DG scheme
				*/
				int ConvDispExIntPattern(std::vector<T>& tripletList, const int offC = 0) {

					/*======================================================*/
					/*			Define Convection Jacobian Block			*/
					/*======================================================*/

					// Convection block [ d RHS_conv / d c ], also depends on upwind entry

					if (_curVelocity >= 0.0) { // forward flow upwind entry -> last node of previous element
						for (unsigned int comp = 0; comp < _nComp; comp++) {
							for (unsigned int i = 0; i < _nNodes; i++) {
								//tripletList.push_back(T(offC + comp * strideColComp() + i * strideColNode(), comp * strideColComp(), 0.0)); // inlet DOFs not included in Jacobian
								for (unsigned int j = 1; j < _nNodes + 1; j++) {
									tripletList.push_back(T(offC + comp * strideColComp() + i * strideColNode(),
										offC + comp * strideColComp() + (j - 1) * strideColNode(),
										0.0));
								}
							}
						}
						for (unsigned int element = 1; element < _nElem; element++) {
							for (unsigned int comp = 0; comp < _nComp; comp++) {
								for (unsigned int i = 0; i < _nNodes; i++) {
									for (unsigned int j = 0; j < _nNodes + 1; j++) {
										// row: jump over inlet DOFs and previous elements, add component offset and go node strides from there for each convection block entry
										// col: jump over inlet DOFs and previous elements, go back one node, add component offset and go node strides from there for each convection block entry
										tripletList.push_back(T(offC + element * strideColElement() + comp * strideColComp() + i * strideColNode(),
											offC + element * strideColElement() - strideColNode() + comp * strideColComp() + j * strideColNode(),
											0.0));
									}
								}
							}
						}
					}
					else { // backward flow upwind entry -> first node of subsequent element
						// special inlet DOF treatment for inlet boundary element (last element)
						for (unsigned int comp = 0; comp < _nComp; comp++) {
							for (unsigned int i = 0; i < _nNodes; i++) {
								// inlet DOFs not included in Jacobian
								for (unsigned int j = 0; j < _nNodes; j++) {
									tripletList.push_back(T(offC + (_nElem - 1) * strideColElement() + comp * strideColComp() + i * strideColNode(),
										offC + (_nElem - 1) * strideColElement() + comp * strideColComp() + j * strideColNode(),
										0.0));
								}
							}
						}
						for (unsigned int element = 0; element < _nElem - 1u; element++) {
							for (unsigned int comp = 0; comp < _nComp; comp++) {
								for (unsigned int i = 0; i < _nNodes; i++) {
									for (unsigned int j = 0; j < _nNodes + 1; j++) {
										// row: jump over inlet DOFs and previous elements, add component offset and go node strides from there for each convection block entry
										// col: jump over inlet DOFs and previous elements, add component offset and go node strides from there for each convection block entry
										tripletList.push_back(T(offC + element * strideColElement() + comp * strideColComp() + i * strideColNode(),
											offC + element * strideColElement() + comp * strideColComp() + j * strideColNode(),
											0.0));
									}
								}
							}
						}
					}

					/*======================================================*/
					/*			Define Dispersion Jacobian Block			*/
					/*======================================================*/

					/* Inner elements */
					if (_nElem >= 5u) {
						// Inner dispersion block [ d RHS_disp / d c ], depends on whole previous and subsequent element plus first entries of subsubsequent elements
						for (unsigned int element = 2; element < _nElem - 2; element++) {
							for (unsigned int comp = 0; comp < _nComp; comp++) {
								for (unsigned int i = 0; i < _nNodes; i++) {
									for (unsigned int j = 0; j < 3 * _nNodes + 2; j++) {
										// row: jump over inlet DOFs and previous elements, add component offset and go node strides from there for each dispersion block entry
										// col: jump over inlet DOFs and previous elements, go back one element and one node, add component offset and go node strides from there for each dispersion block entry
										tripletList.push_back(T(offC + element * strideColElement() + comp * strideColComp() + i * strideColNode(),
											offC + element * strideColElement() - (_nNodes + 1) * strideColNode() + comp * strideColComp() + j * strideColNode(),
											0.0));
									}
								}
							}
						}
					}

					/*		boundary element neighbours		*/

					// left boundary element neighbour
					if (_nElem >= 4u) {
						for (unsigned int comp = 0; comp < _nComp; comp++) {
							for (unsigned int i = 0; i < _nNodes; i++) {
								for (unsigned int j = 1; j < 3 * _nNodes + 2; j++) {
									// row: jump over inlet DOFs and previous element, add component offset and go node strides from there for each dispersion block entry
									// col: jump over inlet DOFs, add component offset and go node strides from there for each dispersion block entry. Also adjust for iterator j (-1)
									tripletList.push_back(T(offC + _nNodes * strideColNode() + comp * strideColComp() + i * strideColNode(),
										offC + comp * strideColComp() + (j - 1) * strideColNode(),
										0.0));
								}
							}
						}
					}
					else if (_nElem == 3u) { // special case: only depends on the two neighbouring elements
						for (unsigned int comp = 0; comp < _nComp; comp++) {
							for (unsigned int i = 0; i < _nNodes; i++) {
								for (unsigned int j = 1; j < 3 * _nNodes + 1; j++) {
									// row: jump over inlet DOFs and previous element, add component offset and go node strides from there for each dispersion block entry
									// col: jump over inlet DOFs, add component offset and go node strides from there for each dispersion block entry. Also adjust for iterator j (-1)
									tripletList.push_back(T(offC + _nNodes * strideColNode() + comp * strideColComp() + i * strideColNode(),
										offC + comp * strideColComp() + (j - 1) * strideColNode(),
										0.0));
								}
							}
						}
					}
					// right boundary element neighbour
					if (_nElem >= 4u) {
						for (unsigned int comp = 0; comp < _nComp; comp++) {
							for (unsigned int i = 0; i < _nNodes; i++) {
								for (unsigned int j = 0; j < 3 * _nNodes + 2 - 1; j++) {
									// row: jump over inlet DOFs and previous elements, add component offset and go node strides from there for each dispersion block entry
									// col: jump over inlet DOFs and previous elements, go back one element and one node, add component offset and go node strides from there for each dispersion block entry.
									tripletList.push_back(T(offC + (_nElem - 2) * strideColElement() + comp * strideColComp() + i * strideColNode(),
										offC + (_nElem - 2) * strideColElement() - (_nNodes + 1) * strideColNode() + comp * strideColComp() + j * strideColNode(),
										0.0));
								}
							}
						}
					}
					/*			boundary elements			*/

					// left boundary element
					unsigned int end = 3u * _nNodes + 2u;
					if (_nElem == 1u) end = 2u * _nNodes + 1u;
					else if (_nElem == 2u) end = 3u * _nNodes + 1u;
					for (unsigned int comp = 0; comp < _nComp; comp++) {
						for (unsigned int i = 0; i < _nNodes; i++) {
							for (unsigned int j = _nNodes + 1; j < end; j++) {
								// row: jump over inlet DOFs, add component offset and go node strides from there for each dispersion block entry
								// col: jump over inlet DOFs, add component offset, adjust for iterator j (-Nnodes-1) and go node strides from there for each dispersion block entry.
								tripletList.push_back(T(offC + comp * strideColComp() + i * strideColNode(),
									offC + comp * strideColComp() + (j - (_nNodes + 1)) * strideColNode(),
									0.0));
							}
						}
					}
					// right boundary element
					if (_nElem >= 3u) {
						for (unsigned int comp = 0; comp < _nComp; comp++) {
							for (unsigned int i = 0; i < _nNodes; i++) {
								for (unsigned int j = 0; j < 2 * _nNodes + 1; j++) {
									// row: jump over inlet DOFs and previous elements, add component offset and go node strides from there for each dispersion block entry
									// col: jump over inlet DOFs and previous elements, go back one element and one node, add component offset and go node strides from there for each dispersion block entry.
									tripletList.push_back(T(offC + (_nElem - 1) * strideColElement() + comp * strideColComp() + i * strideColNode(),
										offC + (_nElem - 1) * strideColElement() - (_nNodes + 1) * strideColNode() + comp * strideColComp() + j * strideColNode(),
										0.0));
								}
							}
						}
					}
					else if (_nElem == 2u) { // special case for _nElem == 2: depends only on left element
						for (unsigned int comp = 0; comp < _nComp; comp++) {
							for (unsigned int i = 0; i < _nNodes; i++) {
								for (unsigned int j = 0; j < 2 * _nNodes; j++) {
									// row: jump over inlet DOFs and previous elements, add component offset and go node strides from there for each dispersion block entry
									// col: jump over inlet DOFs and previous elements, go back one element, add component offset and go node strides from there for each dispersion block entry.
									tripletList.push_back(T(offC + (_nElem - 1) * strideColElement() + comp * strideColComp() + i * strideColNode(),
										offC + (_nElem - 1) * strideColElement() - _nNodes * strideColNode() + comp * strideColComp() + j * strideColNode(),
										0.0));
								}
							}
						}
					}

					return 0;
				}
				/**
				* @brief analytically calculates the convection dispersion jacobian for the collocation DG scheme
				*/
				int calcConvDispCollocationDGSEMJacobian(Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian, Eigen::MatrixXd& jacInlet, const int offC = 0) {

					const int strideColBound = strideColNode() - _nComp;

					/*======================================================*/
					/*			Compute Dispersion Jacobian Block			*/
					/*======================================================*/

					/*		Inner element dispersion blocks		*/

					if (_nElem >= 3u) {
						MatrixXd dispBlock = _DGjacAxDispBlocks[1];
						linalg::BandedEigenSparseRowIterator jacIt(jacobian, offC + strideColElement()); // row iterator starting at second element and component

						for (unsigned int element = 1; element < _nElem - 1; element++) {
							for (unsigned int i = 0; i < dispBlock.rows(); i++, jacIt += strideColBound) {
								for (unsigned int comp = 0; comp < _nComp; comp++, ++jacIt) {
									for (unsigned int j = 0; j < dispBlock.cols(); j++) {
										// pattern is more sparse than a nNodes x 3*nNodes block.
										if ((j >= _nNodes - 1 && j <= 2 * _nNodes) ||
											(i == 0 && j <= 2 * _nNodes) ||
											(i == _nNodes - 1 && j >= _nNodes - 1))
											// row: iterator is at current node i and current component comp
											// col: start at previous element and jump to node j
											jacIt[-strideColElement() + (j - i) * strideColNode()] = dispBlock(i, j) * static_cast<double>(currentDispersion(_curSection)[comp]);
									}
								}
							}
						}
					}

					/*				Boundary element Dispersion blocks			*/

					/* left element */
					MatrixXd dispBlock = _DGjacAxDispBlocks[0];

					if (_nElem != 1u) { // "standard" case
						linalg::BandedEigenSparseRowIterator jacIt(jacobian, offC); // row iterator starting at first element and component

						for (unsigned int i = 0; i < dispBlock.rows(); i++, jacIt += strideColBound) {
							for (unsigned int comp = 0; comp < _nComp; comp++, ++jacIt) {
								for (unsigned int j = _nNodes; j < dispBlock.cols(); j++) {
									// pattern is more sparse than a nNodes x 2*nNodes block.
									if ((j >= _nNodes - 1 && j <= 2 * _nNodes) ||
										(i == 0 && j <= 2 * _nNodes) ||
										(i == _nNodes - 1 && j >= _nNodes - 1))
										// row: iterator is at current node i and current component comp
										// col: jump to node j
										jacIt[((j - _nNodes) - i) * strideColNode()] = dispBlock(i, j) * static_cast<double>(currentDispersion(_curSection)[comp]);
								}
							}
						}
					}
					else { // special case
						linalg::BandedEigenSparseRowIterator jacIt(jacobian, offC); // row iterator starting at first element and component
						for (unsigned int i = 0; i < dispBlock.rows(); i++, jacIt += strideColBound) {
							for (unsigned int comp = 0; comp < _nComp; comp++, ++jacIt) {
								for (unsigned int j = _nNodes; j < _nNodes * 2u; j++) {
									// row: iterator is at current node i and current component comp
									// col: jump to node j
									jacIt[((j - _nNodes) - i) * strideColNode()] = dispBlock(i, j) * static_cast<double>(currentDispersion(_curSection)[comp]);
								}
							}
						}
					}

					/* right element */
					if (_nElem != 1u) { // "standard" case
						dispBlock = _DGjacAxDispBlocks[std::min(_nElem, 3u) - 1];
						linalg::BandedEigenSparseRowIterator jacIt(jacobian, offC + (_nElem - 1) * strideColElement()); // row iterator starting at last element

						for (unsigned int i = 0; i < dispBlock.rows(); i++, jacIt += strideColBound) {
							for (unsigned int comp = 0; comp < _nComp; comp++, ++jacIt) {
								for (unsigned int j = 0; j < 2 * _nNodes; j++) {
									// pattern is more sparse than a nNodes x 2*nNodes block.
									if ((j >= _nNodes - 1 && j <= 2 * _nNodes) ||
										(i == 0 && j <= 2 * _nNodes) ||
										(i == _nNodes - 1 && j >= _nNodes - 1))
										// row: iterator is at current node i and current component comp
										// col: start at previous element and jump to node j
										jacIt[-strideColElement() + (j - i) * strideColNode()] = dispBlock(i, j) * static_cast<double>(currentDispersion(_curSection)[comp]);
								}
							}
						}
					}

					/*======================================================*/
					/*			Compute Convection Jacobian Block			*/
					/*======================================================*/

					// Convection block [ d RHS_conv / d c ], also depends on first entry of previous element
					MatrixXd convBlock = _DGjacAxConvBlock;
					linalg::BandedEigenSparseRowIterator jacIt(jacobian, offC); // row iterator starting at first element and component

					if (_curVelocity > 0.0) { // forward flow upwind convection
						// special inlet DOF treatment for first element (inlet boundary element)
						jacInlet(0, 0) = static_cast<double>(_curVelocity) * convBlock(0, 0); // only first node depends on inlet concentration
						for (unsigned int i = 0; i < convBlock.rows(); i++, jacIt += strideColBound) {
							for (unsigned int comp = 0; comp < _nComp; comp++, ++jacIt) {
								//jacIt[0] = -convBlock(i, 0); // dependency on inlet DOFs is handled in _jacInlet
								for (unsigned int j = 1; j < convBlock.cols(); j++) {
									jacIt[(j - 1 - i) * strideColNode()] += static_cast<double>(_curVelocity) * convBlock(i, j);
								}
							}
						}
						// remaining elements
						for (unsigned int element = 1; element < _nElem; element++) {
							for (unsigned int i = 0; i < convBlock.rows(); i++, jacIt += strideColBound) {
								for (unsigned int comp = 0; comp < _nComp; comp++, ++jacIt) {
									for (unsigned int j = 0; j < convBlock.cols(); j++) {
										// row: iterator is at current element and component
										// col: start at previous elements last node and go to node j.
										jacIt[(j - 1 - i) * strideColNode()] += static_cast<double>(_curVelocity) * convBlock(i, j);
									}
								}
							}
						}
					}
					else { // backward flow upwind convection
						// non-inlet elements
						for (unsigned int element = 0; element < _nElem - 1u; element++) {
							for (unsigned int i = 0; i < convBlock.rows(); i++, jacIt += strideColBound) {
								for (unsigned int comp = 0; comp < _nComp; comp++, ++jacIt) {
									for (unsigned int j = 0; j < convBlock.cols(); j++) {
										// row: iterator is at current element and component
										// col: start at current elements first node and go to node j.
										jacIt[(j - i) * strideColNode()] += static_cast<double>(_curVelocity) * convBlock(i, j);
									}
								}
							}
						}
						// special inlet DOF treatment for last element (inlet boundary element)
						jacInlet(0, 0) = static_cast<double>(_curVelocity) * convBlock(convBlock.rows() - 1, convBlock.cols() - 1); // only last node depends on inlet concentration
						for (unsigned int i = 0; i < convBlock.rows(); i++, jacIt += strideColBound) {
							for (unsigned int comp = 0; comp < _nComp; comp++, ++jacIt) {
								for (unsigned int j = 0; j < convBlock.cols() - 1; j++) {
									jacIt[(j - i) * strideColNode()] += static_cast<double>(_curVelocity) * convBlock(i, j);
								}
							}
						}
					}

					return 0;
				}
				/**
				 * @brief inserts a liquid state block with different factors for components into the system jacobian
				 * @param [in] block (sub)block to be added
				 * @param [in] jac row iterator at first (i.e. upper) entry
				 * @param [in] offCol column to row offset (i.e. start at upper left corner of block)
				 * @param [in] nelements determines how often the block is added (diagonally)
				 * @param [in] Compfactor component dependend factors
				 */
				void insertCompDepLiquidJacBlock(Eigen::MatrixXd block, linalg::BandedEigenSparseRowIterator& jac, int offCol, unsigned int nelements, const active* Compfactor) {

					const int strideColBound = strideColNode() - _nComp;

					for (unsigned int element = 0; element < nelements; element++) {
						for (unsigned int i = 0; i < block.rows(); i++, jac += strideColBound) {
							for (unsigned int comp = 0; comp < _nComp; comp++, ++jac) {
								for (unsigned int j = 0; j < block.cols(); j++) {
									// row: at current node component
									// col: jump to node j
									jac[(j - i) * strideColNode() + offCol] = block(i, j) * static_cast<double>(Compfactor[comp]);
								}
							}
						}
					}
				}
				/**
				 * @brief adds liquid state blocks for all components to the system jacobian
				 * @param [in] block to be added
				 * @param [in] jac row iterator at first (i.e. upper left) entry
				 * @param [in] column to row offset (i.e. start at upper left corner of block)
				 * @param [in] nelements determines how often the block is added (diagonally)
				 */
				void addLiquidJacBlock(Eigen::MatrixXd block, linalg::BandedEigenSparseRowIterator& jac, int offCol, unsigned int nelements) {

					unsigned int strideColBound = strideColNode() - _nComp;

					for (unsigned int element = 0; element < nelements; element++) {
						for (unsigned int i = 0; i < block.rows(); i++, jac += strideColBound) {
							for (unsigned int comp = 0; comp < _nComp; comp++, ++jac) {
								for (unsigned int j = 0; j < block.cols(); j++) {
									// row: at current node component
									// col: jump to node j
									jac[(j - i) * strideColNode() + offCol] += block(i, j);
								}
							}
						}
					}
				}
				/**
				 * @brief analytically calculates the convection dispersion jacobian for the exact integration DG scheme
				 */
				int calcConvDispExIntDGSEMJacobian(Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian, Eigen::MatrixXd& jacInlet, const int offC = 0) {

					/*======================================================*/
					/*			Compute Dispersion Jacobian Block			*/
					/*======================================================*/

					/* Inner elements (exist only if nelements >= 5) */
					if (_nElem >= 5) {
						linalg::BandedEigenSparseRowIterator jacIt(jacobian, offC + strideColElement() * 2); // row iterator starting at third element, first component
						// insert all (nElem - 4) inner element blocks
						insertCompDepLiquidJacBlock(_DGjacAxDispBlocks[2], jacIt, -(strideColElement() + strideColNode()), _nElem - 4u, currentDispersion(_curSection));
					}

					/*	boundary element neighbours (exist only if nelements >= 4)	*/
					if (_nElem >= 4) {
						linalg::BandedEigenSparseRowIterator jacIt(jacobian, offC + strideColElement()); // row iterator starting at second element, first component

						insertCompDepLiquidJacBlock(_DGjacAxDispBlocks[1].block(0, 1, _nNodes, 3 * _nNodes + 1), jacIt, -strideColElement(), 1u, currentDispersion(_curSection));

						jacIt += (_nElem - 4) * strideColElement(); // move iterator to preultimate element (already at third element)
						insertCompDepLiquidJacBlock(_DGjacAxDispBlocks[_nElem > 4 ? 3 : 2].block(0, 0, _nNodes, 3 * _nNodes + 1), jacIt, -(strideColElement() + strideColNode()), 1u, currentDispersion(_curSection));
					}

					/*			boundary elements (exist only if nelements >= 3)			*/
					if (_nElem >= 3) {

						linalg::BandedEigenSparseRowIterator jacIt(jacobian, offC); // row iterator starting at first element, first component

						insertCompDepLiquidJacBlock(_DGjacAxDispBlocks[0].block(0, _nNodes + 1, _nNodes, 2 * _nNodes + 1), jacIt, 0, 1u, currentDispersion(_curSection));

						jacIt += (_nElem - 2) * strideColElement(); // move iterator to last element (already at second element)
						insertCompDepLiquidJacBlock(_DGjacAxDispBlocks[std::min(_nElem, 5u) - 1u].block(0, 0, _nNodes, 2 * _nNodes + 1), jacIt, -(strideColElement() + strideColNode()), 1u, currentDispersion(_curSection));
					}

					/* For special cases nelements = 1, 2, 3, some elements still have to be treated separately*/

					if (_nElem == 1) {
						linalg::BandedEigenSparseRowIterator jacIt(jacobian, offC); // row iterator starting at first element, first component
						// insert the only block
						insertCompDepLiquidJacBlock(_DGjacAxDispBlocks[0].block(0, _nNodes + 1, _nNodes, _nNodes), jacIt, 0, 1u, currentDispersion(_curSection));
					}
					else if (_nElem == 2) {
						linalg::BandedEigenSparseRowIterator jacIt(jacobian, offC); // row iterator starting at first element, first component
						// left Bacobian block
						insertCompDepLiquidJacBlock(_DGjacAxDispBlocks[0].block(0, _nNodes + 1, _nNodes, 2 * _nNodes), jacIt, 0, 1u, currentDispersion(_curSection));
						// right Bacobian block, iterator is already moved to second element
						insertCompDepLiquidJacBlock(_DGjacAxDispBlocks[1].block(0, 1, _nNodes, 2 * _nNodes), jacIt, -strideColElement(), 1u, currentDispersion(_curSection));
					}
					else if (_nElem == 3) {
						linalg::BandedEigenSparseRowIterator jacIt(jacobian, offC + strideColElement()); // row iterator starting at first element, first component
						insertCompDepLiquidJacBlock(_DGjacAxDispBlocks[1].block(0, 1, _nNodes, 3 * _nNodes), jacIt, -strideColElement(), 1u, currentDispersion(_curSection));
					}

					/*======================================================*/
					/*			Compute Convection Jacobian Block			*/
					/*======================================================*/

					linalg::BandedEigenSparseRowIterator jac(jacobian, offC);

					if (_curVelocity > 0.0) { // Forward flow
						// special inlet DOF treatment for inlet (first) element
						jacInlet = static_cast<double>(_curVelocity) * _DGjacAxConvBlock.col(0); // only first element depends on inlet concentration
						addLiquidJacBlock(static_cast<double>(_curVelocity) * _DGjacAxConvBlock.block(0, 1, _nNodes, _nNodes), jac, 0, 1);
						if (_nElem > 1) // iterator already moved to second element
							addLiquidJacBlock(static_cast<double>(_curVelocity) * _DGjacAxConvBlock, jac, -strideColNode(), _nElem - 1);
					}
					else { // Backward flow
						// non-inlet elements first
						if (_nElem > 1)
							addLiquidJacBlock(static_cast<double>(_curVelocity) * _DGjacAxConvBlock, jac, 0, _nElem - 1);
						// special inlet DOF treatment for inlet (last) element. Iterator already moved to last element
						jacInlet = static_cast<double>(_curVelocity) * _DGjacAxConvBlock.col(_DGjacAxConvBlock.cols() - 1); // only last element depends on inlet concentration
						addLiquidJacBlock(static_cast<double>(_curVelocity) * _DGjacAxConvBlock.block(0, 0, _nNodes, _nNodes), jac, 0, 1);
					}

					return 0;
				}
				// todo time derivative pattern
				///**
				//* @brief adds time derivative to the jacobian
				//*/
				//void addTimederJacobian(Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian, double alpha) {

				//	unsigned int offC = 0; // inlet DOFs not included in Jacobian

				//	// =================================================================================================== //
				//	//	 Time derivative Jacobian: d Residual / d y_t   												   //
				//	// =================================================================================================== //

				//	linalg::BandedEigenSparseRowIterator jac(jacobian, offC);

				//	double Beta = (1.0 - static_cast<double>(_totalPorosity)) / static_cast<double>(_totalPorosity);

				//	for (unsigned int point = 0; point < _disc.nPoints; point++) {
				//		for (unsigned int comp = 0; comp < _disc.nComp; comp++) {
				//			// d convDispRHS / d c_t
				//			jac[0] += alpha;

				//			if (_disc.nBound[comp]) { // either one or null; no loop necessary
				//				// d convDispRHS / d q_t
				//				jac[idxr.strideColLiquid() - comp + idxr.offsetBoundComp(comp)] += alpha * Beta;
				//			}
				//			++jac;
				//		}
				//		for (unsigned int comp = 0; comp < _disc.nComp; comp++) {
				//			if (_disc.nBound[comp]) { // either one or null; no loop over bound states necessary
				//				if (!_binding[0]->hasQuasiStationaryReactions()) {
				//					// d isotherm / d q_t
				//					jac[0] += alpha;
				//				}
				//				++jac;
				//			}
				//		}
				//	}
				//}
			};


			/**
			 * @brief Radial convection dispersion transport operator using DG discretization
			 * @details Implements the radial transport equation
			 *
			 * @f[\begin{align}
				\frac{\partial c_i}{\partial t} &= -\frac{u}{\rho} \frac{\partial c_i}{\partial \rho} + \frac{1}{\rho} \frac{\partial}{\partial \rho}\left( D_{\rho,i} \rho \frac{\partial c_i}{\partial \rho} \right) \\
			\end{align} @f]
			 * with Danckwerts boundary conditions
			 * @f[ \begin{align}
			 u c_{\text{in},i}(t) &= u c_i(t,\rho_{in}) - D_{\rho,i} \frac{\partial c_i}{\partial \rho}(t,\rho_{in}) \\
			 \frac{\partial c_i}{\partial \rho}(t,\rho_{out}) &= 0
			 \end{align} @f]
			 *
			 * Key difference from axial: Uses cell-dependent weighted mass matrices M_rho[i] = (delta_rho/2) * M^{(0,1)} + rho_i * M^{(0,0)}
			 *
			 * This class does not store the Jacobian. It only fills existing matrices given to its residual() functions.
			 */
			class RadialConvectionDispersionOperatorBaseDG
			{
			public:

				RadialConvectionDispersionOperatorBaseDG();
				~RadialConvectionDispersionOperatorBaseDG() CADET_NOEXCEPT;

				void setFlowRates(const active& in, const active& out, const active& colPorosity) CADET_NOEXCEPT;

				bool configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, unsigned int nComp, int polynomial_integration_mode, unsigned int nelements, unsigned int polyDeg, unsigned int strideNode);
				bool configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters);
				bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, Eigen::MatrixXd& jacInlet);

				int residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, double* res, Eigen::SparseMatrix<double, Eigen::RowMajor>& jac);
				int residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, active* res, Eigen::SparseMatrix<double, Eigen::RowMajor>& jac);
				int residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, double* res, WithoutParamSensitivity);
				int residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, active* res, WithParamSensitivity);
				int residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, WithParamSensitivity);
				int residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, WithoutParamSensitivity);

				int calcTransportJacobian(Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian, Eigen::MatrixXd& jacInlet, const int bulkOffset = 0);
				typedef Eigen::Triplet<double> T;
				void convDispJacPattern(std::vector<T>& tripletList, const int bulkOffset = 0);
				unsigned int nConvDispEntries(bool pureNNZ = false);
				void multiplyWithDerivativeJacobian(const SimulationTime& simTime, double const* sDot, double* ret) const;
				void addTimeDerivativeToJacobian(double alpha, Eigen::SparseMatrix<double, Eigen::RowMajor>& jacDisc, unsigned int blockOffset = 0);

				// Geometry accessors
				inline const active& columnLength() const CADET_NOEXCEPT { return _colLength; }
				inline const active& innerRadius() const CADET_NOEXCEPT { return _innerRadius; }
				inline const active& outerRadius() const CADET_NOEXCEPT { return _outerRadius; }
				active currentVelocity(double pos) const CADET_NOEXCEPT;
				inline active currentVelocity() const CADET_NOEXCEPT { return currentVelocity(0.0); }
				inline bool forwardFlow() const CADET_NOEXCEPT { return _curVelocity >= 0.0; }

				inline double elemLeftBound(unsigned int idx) const CADET_NOEXCEPT
				{
					return static_cast<double>(_innerRadius) + idx * static_cast<double>(_deltaRho);
				}
				double relativeCoordinate(unsigned int idx) const CADET_NOEXCEPT
				{
					const double rho = static_cast<double>(_innerRadius) +
						(floor(idx / _nNodes) * static_cast<double>(_deltaRho) +
						 0.5 * static_cast<double>(_deltaRho) * (1.0 + _nodes[idx % _nNodes]));
					return (rho - static_cast<double>(_innerRadius)) /
						   (static_cast<double>(_outerRadius) - static_cast<double>(_innerRadius));
				}

				inline const double* LGLnodes() const CADET_NOEXCEPT { return &_nodes[0]; }
				inline const active& currentVelocityCoeff() const CADET_NOEXCEPT { return _curVelocity; }
				inline const active* currentDispersion(const int secIdx) const CADET_NOEXCEPT { return getSectionDependentSlice(_colDispersion, _nComp, secIdx); }
				inline const bool dispersionCompIndep() const CADET_NOEXCEPT { return _dispersionCompIndep; }

				inline unsigned int nComp() const CADET_NOEXCEPT { return _nComp; }
				inline unsigned int nelements() const CADET_NOEXCEPT { return _nElem; }
				inline unsigned int nNodes() const CADET_NOEXCEPT { return _nNodes; }
				inline unsigned int nPoints() const CADET_NOEXCEPT { return _nPoints; }
				inline bool exactInt() const CADET_NOEXCEPT { return true; }  // Radial DG always uses exact integration

				// Indexer functionality
				inline int strideColElement() const CADET_NOEXCEPT { return static_cast<int>(_strideElem); }
				inline int strideColNode() const CADET_NOEXCEPT { return static_cast<int>(_strideNode); }
				inline int strideColComp() const CADET_NOEXCEPT { return 1; }
				inline int offsetC() const CADET_NOEXCEPT { return _nComp; }

				unsigned int jacobianLowerBandwidth() const CADET_NOEXCEPT;
				unsigned int jacobianUpperBandwidth() const CADET_NOEXCEPT;
				double inletJacobianFactor() const CADET_NOEXCEPT;
				unsigned int nJacEntries(bool pureNNZ);

				// Radial DG always uses exact integration (4 * nNodes stencil)
				int requiredADdirs() const CADET_NOEXCEPT { return 4 * _nNodes * strideColNode() + 1; }

				// Accessors for mass matrices (needed for rLRMP film diffusion)
				inline const Eigen::MatrixXd& invMRho(unsigned int cell) const CADET_NOEXCEPT { return _invMM_rho[cell]; }
				inline const Eigen::MatrixXd& M00() const CADET_NOEXCEPT { return _M00; }
				inline const Eigen::MatrixXd& M01() const CADET_NOEXCEPT { return _M01; }
				inline double deltaRho() const CADET_NOEXCEPT { return static_cast<double>(_deltaRho); }

				/**
				 * @brief Computes the weighted mass matrix M_rho for a given cell
				 * @detail M_rho = (delta_rho/2) * M01 + rho_left * M00
				 */
				inline Eigen::MatrixXd MRho(unsigned int cell) const
				{
					const double rho_left = elemLeftBound(cell);
					return (static_cast<double>(_deltaRho) / 2.0) * _M01 + rho_left * _M00;
				}

				bool setParameter(const ParameterId& pId, double value);
				bool setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue);
				bool setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& id, double value);

			protected:

				// Discretization parameters
				unsigned int _nComp;      //!< Number of components
				unsigned int _polyDeg;    //!< DG polynomial degree
				unsigned int _nElem;      //!< Number of radial elements
				unsigned int _nNodes;     //!< Nodes per element
				unsigned int _nPoints;    //!< Total radial discrete points

				unsigned int _strideNode; //!< Stride between nodes
				unsigned int _strideElem; //!< Stride between elements

				// Radial geometry
				active _colLength;        //!< Column height L
				active _innerRadius;      //!< Inner radius rho_inner
				active _outerRadius;      //!< Outer radius rho_outer
				active _deltaRho;         //!< Element spacing (uniform grid)

				// Standard DG matrices (in reference space)
				Eigen::VectorXd _nodes;       //!< LGL nodes in [-1, 1]
				Eigen::MatrixXd _polyDerM;    //!< Polynomial derivative matrix D
				Eigen::VectorXd _invWeights;  //!< Inverse LGL quadrature weights (collocation)
				Eigen::MatrixXd _invMM;       //!< Inverse standard mass matrix M^{(0,0)}^{-1}
				Eigen::MatrixXd _M00;         //!< Standard mass matrix M^{(0,0)}
				Eigen::MatrixXd _M01;         //!< Weighted mass matrix M^{(0,1)}

				// Cell-dependent matrices for radial geometry
				std::vector<Eigen::MatrixXd> _invMM_rho;    //!< Per-cell inverse weighted mass matrix M_rho^{-1}[cell]
				std::vector<Eigen::MatrixXd> _S_g;          //!< Per-cell stiffness matrix D^T * M_rho [cell]
				std::vector<double> _rhoNodeCoords;         //!< Physical radial coordinates of all DG nodes
				std::vector<double> _rhoCellBounds;         //!< Cell boundary positions (rho at interfaces)

				// Per-cell Jacobian blocks
				std::vector<Eigen::MatrixXd> _DGjacRadDispBlocks;  //!< Cell-dependent dispersion Jacobian blocks
				std::vector<Eigen::MatrixXd> _DGjacRadConvBlocks;  //!< Cell-dependent convection Jacobian blocks

				// Auxiliary state for dispersion
				active* _auxState;        //!< Auxiliary variable g = d c / d rho
				active* _subsState;       //!< Substitute variable
				Eigen::Vector<active, Eigen::Dynamic> _surfaceFlux;  //!< Interface flux values
				Eigen::Vector<active, 4> _boundary;  //!< Boundary values

				// Section dependent parameters
				std::vector<active> _colDispersion;  //!< Radial dispersion D_rho
				std::vector<active> _velocity;       //!< Velocity coefficient (section dependent)
				active _curVelocity;                 //!< Current velocity coefficient u = Q/(2*pi*L)
				int _dir;                            //!< Flow direction

				int _curSection;
				bool _newStaticJac;

				bool _dispersionCompIndep;
				IParameterParameterDependence* _dispersionDep;

				// Variable dispersion support
				bool _variableDispersion;                           //!< Flag indicating spatially-varying dispersion
				bool _overintegrate;                                //!< Use 3p/2 quadrature for nonlinear coefficients
				std::vector<Eigen::VectorXd> _dispersionAtNodes;    //!< Per-cell dispersion values at DG nodes [cell][node]
				std::vector<double> _dispersionAtInterfaces;        //!< Dispersion values at cell interfaces [nElem+1]

				/* ===================================================================================
				 *  Helper functions for radial DG
				 * =================================================================================== */

				/**
				 * @brief Computes cell-dependent weighted mass matrices and stiffness matrices
				 * @param [in] dispersion base dispersion coefficient (used when not variable)
				 */
				void computeCellDependentMatrices(double dispersion = 1.0);

				/**
				 * @brief Updates dispersion values at nodes and interfaces for variable D(rho)
				 * @param [in] model reference to the model (needed for parameter dependence evaluation)
				 * @param [in] secIdx current section index
				 * @param [in] comp component index
				 */
				void updateDispersionValues(const IModel& model, unsigned int secIdx, unsigned int comp);

				/**
				 * @brief Recomputes S_g matrices for variable dispersion
				 * @detail Uses quadrature integration with dispersion values at DG nodes
				 */
				void recomputeDispersionMatrices();

				/**
				 * @brief Computes radial node coordinates in physical space
				 */
				void computeRadialNodeCoordinates();

				/**
				 * @brief Calculates the convection part of the radial DG Jacobian for a specific cell
				 * @param [in] cellIdx cell index
				 */
				Eigen::MatrixXd DGjacobianConvBlockRadial(unsigned int cellIdx);

				/**
				 * @brief Calculates the dispersion part of the radial DG Jacobian for a specific cell
				 * @param [in] cellIdx cell index
				 */
				Eigen::MatrixXd DGjacobianDispBlockRadial(unsigned int cellIdx);

				/**
				 * @brief Calculates the auxiliary G block for a specific cell
				 * @param [in] cellIdx cell index
				 */
				Eigen::MatrixXd getGBlockRadial(unsigned int cellIdx);

				/**
				 * @brief Calculates the radial DG Jacobian (exact integration only)
				 */
				void calcConvDispRadialDGSEMJacobian(Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian, Eigen::MatrixXd& jacInlet, const int offC);

				/* ===================================================================================
				 *  Radial DG residual helper functions
				 * =================================================================================== */

				/**
				 * @brief Computes volume integral for radial DG with cell-dependent stiffness
				 */
				template<typename StateType, typename ResidualType>
				void volumeIntegralRadial(unsigned int cellIdx,
					Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state,
					Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>>& stateDer);

				/**
				 * @brief Computes surface integral for radial DG with rho-weighting
				 */
				template<typename StateType, typename ResidualType>
				void surfaceIntegralRadial(unsigned int cellIdx,
					const StateType* state, ResidualType* stateDer,
					unsigned int strideNode_state, unsigned int strideElem_state,
					unsigned int strideNode_stateDer, unsigned int strideElem_stateDer);

				/**
				 * @brief Computes numerical flux for auxiliary equation (central flux)
				 */
				template<typename StateType>
				void InterfaceFluxAuxiliaryRadial(const StateType* C, unsigned int strideNode, unsigned int strideElem);

				/**
				 * @brief Computes numerical flux for main equation (convection + dispersion)
				 */
				template<typename StateType, typename ParamType>
				void InterfaceFluxRadial(const StateType* C, ParamType dispersion);

				/**
				 * @brief Calculates boundary values for radial Danckwerts conditions
				 */
				template<typename StateType>
				void calcBoundaryValuesRadial() {
					// g_l left boundary (inlet for forward flow)
					_boundary[2] = -reinterpret_cast<StateType*>(_auxState)[0];
					// g_r right boundary (outlet for forward flow)
					_boundary[3] = -reinterpret_cast<StateType*>(_auxState)[_nPoints - 1];
				}

				/* ===================================================================================
				 *  Radial DG residual template implementations
				 * =================================================================================== */

				/**
				 * @brief Computes volume integral for auxiliary equation: g -= D * c
				 */
				template<typename StateType, typename ResidualType>
				void volumeIntegralAuxRadial(
					Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state,
					Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>>& stateDer) {

					for (unsigned int elem = 0; elem < _nElem; elem++) {
						if constexpr (std::is_same_v<StateType, double>) {
							if constexpr (std::is_same_v<ResidualType, double>) {
								stateDer.segment(elem * _nNodes, _nNodes) -= _polyDerM * state.segment(elem * _nNodes, _nNodes);
							}
							else {
								stateDer.segment(elem * _nNodes, _nNodes) -= (_polyDerM * state.segment(elem * _nNodes, _nNodes)).template cast<ResidualType>();
							}
						}
						else {
							stateDer.segment(elem * _nNodes, _nNodes) -= _polyDerM.template cast<active>() * state.segment(elem * _nNodes, _nNodes);
						}
					}
				}

				/**
				 * @brief Computes volume integral for main equation with radial stiffness
				 * @detail Uses cell-dependent S_g[cell] = D^T * M_rho for dispersion term
				 */
				template<typename StateType, typename ResidualType, typename ParamType>
				void volumeIntegralMainRadial(const StateType* c, ResidualType* res, const StateType* g, ParamType d_rad) {
					const double deltaRho = static_cast<double>(_deltaRho);

					for (unsigned int elem = 0; elem < _nElem; elem++) {
						// Convection term: M_rho^{-1} * D^T * M^{(0,0)} * u * c
						// Dispersion term: -M_rho^{-1} * S_g * g

						for (unsigned int node = 0; node < _nNodes; node++) {
							ResidualType convSum = 0.0;
							ResidualType dispSum = 0.0;

							for (unsigned int j = 0; j < _nNodes; j++) {
								// Convection: D^T * M00 * c
								double DtM00_entry = 0.0;
								for (unsigned int k = 0; k < _nNodes; k++) {
									DtM00_entry += _polyDerM(k, node) * _M00(k, j);
								}
								convSum += static_cast<ResidualType>(_invMM_rho[elem](node, j) * DtM00_entry * _curVelocity * c[elem * _strideElem + j * _strideNode]);

								// Dispersion: S_g * g (S_g = D^T * M_rho, already computed)
								dispSum += static_cast<ResidualType>(_invMM_rho[elem](node, j) * _S_g[elem](j, j) * d_rad * g[elem * _nNodes + j]);
							}

							res[elem * _strideElem + node * _strideNode] += static_cast<ResidualType>((2.0 / deltaRho) * (convSum - dispSum));
						}
					}
				}

				/**
				 * @brief Computes surface integral for auxiliary equation with radial geometry
				 * @detail Exact integration: g -= M^{-1} * B * [c - c*]
				 */
				template<typename StateType, typename ResidualType>
				void surfaceIntegralAuxRadial(const StateType* state, ResidualType* stateDer,
					unsigned int strideNode_state, unsigned int strideElem_state,
					unsigned int strideNode_stateDer, unsigned int strideElem_stateDer) {

					for (unsigned int elem = 0; elem < _nElem; elem++) {
						// M^{-1} B [state - state*]
						for (unsigned int node = 0; node < _nNodes; node++) {
							stateDer[elem * strideElem_stateDer + node * strideNode_stateDer]
								-= static_cast<ResidualType>(_invMM(node, 0) * (state[elem * strideElem_state] - _surfaceFlux[elem])
									- _invMM(node, _polyDeg) * (state[elem * strideElem_state + _polyDeg * strideNode_state] - _surfaceFlux[elem + 1]));
						}
					}
				}

				/**
				 * @brief Computes surface integral for main equation with radial weighting
				 * @detail Exact integration: Dc -= (2/Δρ) * M_ρ^{-1} * B * h*
				 *         Note: This version is not used - surfaceIntegralMainRadialImpl is used instead
				 */
				template<typename StateType, typename ResidualType, typename ParamType>
				void surfaceIntegralMainRadial(ResidualType* res, ParamType d_rad) {
					const double deltaRho = static_cast<double>(_deltaRho);

					for (unsigned int elem = 0; elem < _nElem; elem++) {
						for (unsigned int node = 0; node < _nNodes; node++) {
							// Left face: -M_rho^{-1}[:, 0] * h*_left
							res[elem * _strideElem + node * _strideNode]
								-= static_cast<ResidualType>((2.0 / deltaRho) * _invMM_rho[elem](node, 0) * _surfaceFlux[elem]);

							// Right face: +M_rho^{-1}[:, N] * h*_right
							res[elem * _strideElem + node * _strideNode]
								+= static_cast<ResidualType>((2.0 / deltaRho) * _invMM_rho[elem](node, _nNodes - 1) * _surfaceFlux[elem + 1]);
						}
					}
				}

				/**
				 * @brief Computes auxiliary flux c* (central average)
				 */
				template<typename StateType>
				void InterfaceFluxAuxiliaryRadialImpl(const StateType* C, unsigned int strideNode, unsigned int strideElem) {
					// c* = 0.5 * (c_left + c_right) at interior interfaces
					for (unsigned int elem = 1; elem < _nElem; elem++) {
						_surfaceFlux[elem] = 0.5 * (C[elem * strideElem - strideNode] + C[elem * strideElem]);
					}

					// Boundary fluxes
					_surfaceFlux[0] = C[0];  // Left boundary: just use interior value
					_surfaceFlux[_nElem] = C[_nElem * strideElem - strideNode];  // Right boundary
				}

				/**
				 * @brief Computes numerical fluxes c* and g* for radial DG
				 * @detail Following CADET-Julia pattern:
				 *   - c*: upwind convection flux
				 *   - g*: central average for interior, Danckwerts for inlet, zero for outlet
				 *   - u = v * rho_inlet for forward flow (velocity × boundary radius)
				 * @param [in] C concentration state
				 * @param [in] d_rad base dispersion coefficient (used if not variable dispersion)
				 * @param [out] c_star concentration flux at interfaces
				 * @param [out] g_star gradient flux at interfaces
				 * @param [out] d_rad_i dispersion at interfaces
				 */
				template<typename StateType, typename ParamType>
				void computeNumericalFluxesRadial(const StateType* C, ParamType d_rad,
					std::vector<StateType>& c_star, std::vector<StateType>& g_star, std::vector<double>& d_rad_i) {

					StateType* g = reinterpret_cast<StateType*>(_auxState);
					const double deltaRho = static_cast<double>(_deltaRho);

					unsigned int strideNode = _strideNode;
					unsigned int strideElem = _nNodes * strideNode;
					unsigned int strideNode_g = 1u;
					unsigned int strideElem_g = _nNodes * strideNode_g;

					// Velocity coefficient (NOT multiplied by ρ - matches FV implementation)
					// _curVelocity = Q/(2πL), and actual velocity v(ρ) = _curVelocity/ρ
					// The flux is u_coeff * c (constant), not ρ*v*c
					const double u = static_cast<double>(_curVelocity);

					// Dispersion at interfaces: use pre-computed values if variable, otherwise constant
					if (_variableDispersion) {
						for (unsigned int i = 0; i <= _nElem; i++) {
							d_rad_i[i] = _dispersionAtInterfaces[i];
						}
					} else {
						for (unsigned int i = 0; i <= _nElem; i++) {
							d_rad_i[i] = static_cast<double>(d_rad);
						}
					}

					// Interior interfaces
					for (unsigned int elem = 1; elem < _nElem; elem++) {
						// c*: upwind flux
						if (_curVelocity >= 0.0) {
							c_star[elem] = C[elem * strideElem - strideNode];  // from left
						} else {
							c_star[elem] = C[elem * strideElem];  // from right
						}
						// g*: central average
						g_star[elem] = 0.5 * (g[elem * strideElem_g - strideNode_g] + g[elem * strideElem_g]);
					}

					// Boundary interfaces
					// Note: _boundary[0] is of type 'active', so we use static_cast<StateType> directly
					// to preserve AD tracking when StateType is 'active'
					// g is now in REFERENCE coordinates (dc/dξ), not physical (dc/dρ)
					// Danckwerts BC: u*c_in = u*c - ρ*D*(2/Δρ)*g_ref
					// => g_ref = u*(c - c_in)*Δρ / (2*ρ*D)
					if (_curVelocity >= 0.0) {
						// Inlet at left (rho_inner)
						c_star[0] = static_cast<StateType>(_boundary[0]);  // inlet concentration
						const double rho_inlet = _rhoCellBounds[0];
						if (std::abs(d_rad_i[0]) > 1e-30 && std::abs(rho_inlet) > 1e-30) {
							g_star[0] = static_cast<StateType>(u * deltaRho / (2.0 * rho_inlet * d_rad_i[0]))
							          * (C[0] - static_cast<StateType>(_boundary[0]));
						} else {
							g_star[0] = StateType(0.0);
						}

						// Outlet at right (rho_outer) - zero flux BC (g* = 0)
						c_star[_nElem] = C[_nElem * strideElem - strideNode];
						g_star[_nElem] = StateType(0.0);
					} else {
						// Inlet at right for backward flow
						c_star[_nElem] = static_cast<StateType>(_boundary[0]);
						const double rho_inlet = _rhoCellBounds[_nElem];
						if (std::abs(d_rad_i[_nElem]) > 1e-30 && std::abs(rho_inlet) > 1e-30) {
							g_star[_nElem] = static_cast<StateType>(u * deltaRho / (2.0 * rho_inlet * d_rad_i[_nElem]))
							               * (C[_nElem * strideElem - strideNode] - static_cast<StateType>(_boundary[0]));
						} else {
							g_star[_nElem] = StateType(0.0);
						}

						// Outlet at left for backward flow - zero flux BC
						c_star[0] = C[0];
						g_star[0] = StateType(0.0);
					}
				}

				/**
				 * @brief Computes surface integral for main equation (CADET-Julia pattern)
				 * @detail Exact integration: Dc -= (2/Δρ) * M_ρ^{-1} * [
				 *           [:, 0] * (-u * c*[left] + rho[left] * D[left] * g*[left]) +
				 *           [:, N] * (u * c*[right] - rho[right] * D[right] * g*[right])
				 *         ]
				 * @param [in,out] res pointer to residual for current component
				 * @param [in] strideNode stride between nodes in res
				 * @param [in] strideElem stride between elements in res
				 */
				template<typename StateType, typename ResidualType>
				void surfaceIntegralMainRadialImpl(ResidualType* res,
					const std::vector<StateType>& c_star, const std::vector<StateType>& g_star,
					const std::vector<double>& d_rad_i, double u,
					unsigned int strideNode, unsigned int strideElem) {

					const double deltaRho = static_cast<double>(_deltaRho);

					for (unsigned int elem = 0; elem < _nElem; elem++) {
						double rho_left = _rhoCellBounds[elem];
						double rho_right = _rhoCellBounds[elem + 1];

						// Flux F = u*c - ρ*D*(2/Δρ)*g_ref  (g_star is in reference coords)
						// Left face contribution: -F* = -u*c* + ρ*D*(2/Δρ)*g*
						ResidualType left_flux = static_cast<ResidualType>(
							-u * c_star[elem] + rho_left * d_rad_i[elem] * (2.0 / deltaRho) * g_star[elem]
						);

						// Right face contribution: F* = u*c* - ρ*D*(2/Δρ)*g*
						ResidualType right_flux = static_cast<ResidualType>(
							u * c_star[elem + 1] - rho_right * d_rad_i[elem + 1] * (2.0 / deltaRho) * g_star[elem + 1]
						);

						// Surface integral: res += M_ρ^{-1} * B * F* (strong form correction)
						for (unsigned int node = 0; node < _nNodes; node++) {
							res[elem * strideElem + node * strideNode]
								+= static_cast<ResidualType>(
									_invMM_rho[elem](node, 0) * left_flux +
									_invMM_rho[elem](node, _nNodes - 1) * right_flux
								);
						}
					}
				}

				/**
				 * @brief Computes volume integral for main equation (CADET-Julia pattern)
				 * @detail Dc += (2/Δρ) * M_ρ^{-1} * (D^T * M^{(0,0)} * u * c - S_g * g)
				 *         For constant dispersion: S_g = D^T * M_ρ, multiply by d_rad
				 *         For variable dispersion: S_g already includes D(ρ), don't multiply by d_rad
				 * @param [in] c pointer to concentration for current component
				 * @param [in,out] res pointer to residual for current component
				 * @param [in] g auxiliary variable (single component, stride 1)
				 * @param [in] strideNode_c stride between nodes in c
				 * @param [in] strideElem_c stride between elements in c
				 * @param [in] strideNode_res stride between nodes in res
				 * @param [in] strideElem_res stride between elements in res
				 */
				template<typename StateType, typename ResidualType>
				void volumeIntegralMainRadialImpl(const StateType* c, ResidualType* res, const StateType* g, double u, double d_rad,
					unsigned int strideNode_c, unsigned int strideElem_c,
					unsigned int strideNode_res, unsigned int strideElem_res) {
					const double deltaRho = static_cast<double>(_deltaRho);

					// For variable dispersion, d_rad is already baked into S_g
					const double dispFactor = _variableDispersion ? 1.0 : d_rad;

					for (unsigned int elem = 0; elem < _nElem; elem++) {
						// Precompute D^T * M^{(0,0)} * c for this element
						Eigen::VectorXd DtM00_c = Eigen::VectorXd::Zero(_nNodes);
						for (unsigned int i = 0; i < _nNodes; i++) {
							for (unsigned int j = 0; j < _nNodes; j++) {
								double DtM00_ij = 0.0;
								for (unsigned int k = 0; k < _nNodes; k++) {
									DtM00_ij += _polyDerM(k, i) * _M00(k, j);
								}
								DtM00_c(i) += DtM00_ij * static_cast<double>(c[elem * strideElem_c + j * strideNode_c]);
							}
						}

						// Precompute S_g * g for this element (g has stride 1)
						Eigen::VectorXd Sg_g = Eigen::VectorXd::Zero(_nNodes);
						for (unsigned int i = 0; i < _nNodes; i++) {
							for (unsigned int j = 0; j < _nNodes; j++) {
								Sg_g(i) += _S_g[elem](i, j) * static_cast<double>(g[elem * _nNodes + j]);
							}
						}

						// Apply M_ρ^{-1} and add to residual
						// NOTE: No (2/Δρ) factor here - coordinate transforms cancel in weak form
						// See: ∫(dL/dρ)*F*dρ = ∫(2/Δρ)*(dL/dξ)*F*(Δρ/2)*dξ = ∫(dL/dξ)*F*dξ
						for (unsigned int node = 0; node < _nNodes; node++) {
							double conv_term = 0.0;
							double disp_term = 0.0;

							for (unsigned int j = 0; j < _nNodes; j++) {
								conv_term += _invMM_rho[elem](node, j) * DtM00_c(j);
								disp_term += _invMM_rho[elem](node, j) * Sg_g(j);
							}

							// Sign: res = dc/dt - M_ρ^{-1} * D * F for strong form DG
							res[elem * strideElem_res + node * strideNode_res]
								-= static_cast<ResidualType>(u * conv_term - dispFactor * disp_term);
						}
					}
				}

				/* ===================================================================================
				 *  Main radial DG residualImpl template implementation
				 * =================================================================================== */

				/**
				 * @brief Main residual implementation for radial DG
				 * @detail Following CADET-Julia radialresidualImpl! pattern:
				 *   1. Initialize residual with time derivative dc/dt
				 *   2. For each component:
				 *      a. Compute auxiliary variable g = dc/drho via volume + surface integrals
				 *      b. Scale g by (2/Δρ)
				 *      c. Compute numerical fluxes c*, g*
				 *      d. Compute main equation Dc via surface + volume integrals
				 */
				template <typename StateType, typename ResidualType, typename ParamType, typename RowIteratorType, bool wantJac>
				int residualImpl(const IModel& model, double t, unsigned int secIdx,
					StateType const* y, double const* yDot, ResidualType* res, RowIteratorType jacBegin) {

					const double deltaRho = static_cast<double>(_deltaRho);

					// Working strides
					unsigned int strideNode = _strideNode;
					unsigned int strideElem = _nNodes * strideNode;
					unsigned int strideNode_g = 1u;
					unsigned int strideElem_g = _nNodes;

					// Loop over all components (like axial DG)
					for (unsigned int comp = 0; comp < _nComp; comp++) {

						// Get dispersion coefficient for this component (section-dependent)
						ParamType d_rad = static_cast<ParamType>(getSectionDependentSlice(_colDispersion, _nComp, secIdx)[comp]);

						// Pointer to bulk DOFs for this component (skip inlet DOFs)
						StateType const* yBulk = y + offsetC() + comp;

						// Pointer to bulk residual for this component
						ResidualType* resBulk = res + offsetC() + comp;

						// ============================================================
						// Step 0: Initialize residual with time derivative dc/dt
						// ============================================================
						if (yDot) {
							for (unsigned int i = 0; i < _nPoints; i++) {
								resBulk[i * strideNode] = static_cast<ResidualType>(yDot[offsetC() + comp + i * strideNode]);
							}
						} else {
							for (unsigned int i = 0; i < _nPoints; i++) {
								resBulk[i * strideNode] = ResidualType(0.0);
							}
						}

						// Set inlet boundary concentration for this component
						if constexpr (std::is_same_v<StateType, active>) {
							_boundary[0] = y[comp];  // AD tracking preserved
						} else {
							_boundary[0] = static_cast<active>(y[comp]);
						}

						// Update variable dispersion values if needed
						if (_variableDispersion && _newStaticJac) {
							auto* mutableThis = const_cast<RadialConvectionDispersionOperatorBaseDG*>(this);
							mutableThis->updateDispersionValues(model, secIdx, comp);
							mutableThis->recomputeDispersionMatrices();
						}

						// Auxiliary state g for this component
						StateType* g = reinterpret_cast<StateType*>(_auxState);

						// ============================================================
						// Step 1: Compute auxiliary variable g = dc/drho
						// ============================================================

						// Zero out g
						for (unsigned int i = 0; i < _nPoints; i++) {
							g[i] = StateType(0.0);
						}

						// Create Eigen maps for volume integral
						Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> yMap(
							yBulk, _nPoints, InnerStride<Dynamic>(strideNode));
						Eigen::Map<Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> gMap(
							g, _nPoints, InnerStride<Dynamic>(strideNode_g));

						// Volume integral: g -= D * c
						volumeIntegralAuxRadial<StateType, StateType>(yMap, gMap);

						// Compute auxiliary flux c* (central average for interior)
						InterfaceFluxAuxiliaryRadialImpl<StateType>(yBulk, strideNode, strideElem);

						// Surface integral for g
						surfaceIntegralAuxRadial<StateType, StateType>(yBulk, g, strideNode, strideElem, strideNode_g, strideElem_g);

						// Note: g is in reference coordinates (dc/dξ)
						// The (2/Δρ) scaling for physical gradient is handled implicitly
						// through the S_g matrix which includes the coordinate transform

						// ============================================================
						// Step 2: Compute numerical fluxes c* and g*
						// ============================================================

						std::vector<StateType> c_star(_nElem + 1);
						std::vector<StateType> g_star(_nElem + 1);
						std::vector<double> d_rad_i(_nElem + 1);

						computeNumericalFluxesRadial<StateType, ParamType>(yBulk, d_rad, c_star, g_star, d_rad_i);

						// Velocity coefficient (NOT multiplied by ρ - matches FV implementation)
						const double u = static_cast<double>(_curVelocity);

						// ============================================================
						// Step 3: Compute main equation residual (convection + dispersion)
						// ============================================================

						// Surface integral: res -= (2/Δρ) * M_ρ^{-1} * [flux terms]
						surfaceIntegralMainRadialImpl<StateType, ResidualType>(resBulk, c_star, g_star, d_rad_i, u,
							strideNode, strideElem);

						// Volume integral: res += (2/Δρ) * M_ρ^{-1} * (convection - dispersion)
						volumeIntegralMainRadialImpl<StateType, ResidualType>(yBulk, resBulk, g, u, static_cast<double>(d_rad),
							strideNode, strideElem, strideNode, strideElem);

					} // end component loop

					return 0;
				}
			};

			/**
			 * @brief Convection dispersion transport operator
			 * @details This class wraps AxialConvectionDispersionOperatorBaseDG and provides all the functionality it does.
			 * In addition, the Jacobian is stored and corresponding functions are provided
			 * (assembly, factorization, solution, retrieval).
			 *
			 * This class assumes that the first element is offset by the number of components (inlet DOFs) in the global state vector.
			 */
			template <typename BaseOperator>
			class ConvectionDispersionOperatorDG
			{
			public:

				ConvectionDispersionOperatorDG();
				~ConvectionDispersionOperatorDG() CADET_NOEXCEPT;

				int requiredADdirs() const CADET_NOEXCEPT;

				void setFlowRates(const active& in, const active& out, const active& colPorosity) CADET_NOEXCEPT;

				bool configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, unsigned int nComp, int polynomial_integration_mode, unsigned int nelements, unsigned int polyDeg, unsigned int strideNode);
				bool configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters);
				bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const AdJacobianParams& adJac);

				int residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, double* res, Eigen::SparseMatrix<double, Eigen::RowMajor>& jac);
				int residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, active* res, Eigen::SparseMatrix<double, Eigen::RowMajor>& jac);
				int residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, double* res, bool wantJac, WithoutParamSensitivity);
				int residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, bool wantJac, WithParamSensitivity);
				int residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, bool wantJac, WithoutParamSensitivity);
				int residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, active* res, bool wantJac, WithParamSensitivity);

				void prepareADvectors(const AdJacobianParams& adJac) const;
				void extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset);

				bool solveTimeDerivativeSystem(const SimulationTime& simTime, double* const rhs);
				void multiplyWithDerivativeJacobian(const SimulationTime& simTime, double const* sDot, double* ret) const;

				bool assembleAndFactorizeDiscretizedJacobian(double alpha);
				bool solveDiscretizedJacobian(double* rhs) const;

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
				double checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const;
#endif

				inline const active& columnLength() const CADET_NOEXCEPT { return _baseOp.columnLength(); }
				inline active currentVelocity(double pos) const CADET_NOEXCEPT { return _baseOp.currentVelocity(pos); }
				inline bool forwardFlow() const CADET_NOEXCEPT { return _baseOp.forwardFlow(); }

				inline double elemLeftBound(unsigned int idx) const CADET_NOEXCEPT { return _baseOp.elemLeftBound(idx); }
				inline double relativeCoordinate(unsigned int idx) const CADET_NOEXCEPT { return _baseOp.relativeCoordinate(idx); }

				inline Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian() CADET_NOEXCEPT { return _jacC; }
				inline const Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian() const CADET_NOEXCEPT { return _jacC; }

				inline Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobianDisc() CADET_NOEXCEPT { return _jacCdisc; }
				inline const Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobianDisc() const CADET_NOEXCEPT { return _jacCdisc; }

				inline bool setParameter(const ParameterId& pId, double value)
				{
					return _baseOp.setParameter(pId, value);
				}
				inline bool setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue)
				{
					return _baseOp.setSensitiveParameter(sensParams, pId, adDirection, adValue);
				}
				inline bool setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& id, double value)
				{
					return _baseOp.setSensitiveParameterValue(sensParams, id, value);
				}

			protected:

				void addTimeDerivativeToJacobian(double alpha);
				void assembleDiscretizedJacobian(double alpha);

				BaseOperator _baseOp;

				Eigen::MatrixXd _jacInlet; //!< inlet Jacobian
				Eigen::SparseMatrix<double, Eigen::RowMajor> _jacC; //!< Jacobian
				Eigen::SparseMatrix<double, Eigen::RowMajor> _jacCdisc; //!< Jacobian with time derivatives from BDF method

				// Indexer functionality
				// Offsets
				inline int offsetC() const CADET_NOEXCEPT { return _baseOp.nComp(); }

			};

			extern template class ConvectionDispersionOperatorDG<AxialConvectionDispersionOperatorBaseDG>;
			extern template class ConvectionDispersionOperatorDG<RadialConvectionDispersionOperatorBaseDG>;

			typedef ConvectionDispersionOperatorDG<AxialConvectionDispersionOperatorBaseDG> AxialConvectionDispersionOperatorDG;
			typedef ConvectionDispersionOperatorDG<RadialConvectionDispersionOperatorBaseDG> RadialConvectionDispersionOperatorDG;

		} // namespace parts
	} // namespace model
} // namespace cadet

#endif  // LIBCADET_CONVECTIONDISPERSIONOPERATORDG_HPP_
