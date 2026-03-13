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
 * Defines a discrete calculus toolbox for nodal spectral elements used in Discontinuous Galerkin Spectral Element Methods (DGSEM).
 * @details See @cite Kopriva2009
 */

#ifndef LIBCADET_DGTOOLBOX_HPP_
#define LIBCADET_DGTOOLBOX_HPP_

#include <Eigen/Dense>

namespace cadet
{

namespace model
{

namespace parts
{

namespace dgtoolbox
{
/**
 * @brief computes the node in physical space corresponding to a reference node
 * @param [in] delta element spacing
 * @param [in] elemIdx element index starting at 0
 * @param [in] xi reference node
 */
template <typename ParamType>
ParamType mapRefToPhys(const std::vector<ParamType> delta, const unsigned int elemIdx, const double xi)
{
	//return std::accumulate(delta.begin(), delta.begin() + elemIdx, delta[elemIdx] / 2.0 * (xi + 1.0));

	ParamType map = 0.0;
	for (unsigned int i = 0; i < elemIdx; i++)
		map += delta[i];
	return map + delta[elemIdx] / 2.0 * (xi + 1.0);
}
/**
 * @brief computes the node in referece space corresponding to a physical node
 * @param [in] deltaX element spacing
 * @param [in] elemIdx element index starting at 0
 * @param [in] x physical node
 */
template <typename ParamType>
ParamType mapPhysToRef(const std::vector<ParamType> deltaX, const unsigned int elemIdx, const double x)
{
	//return (x - std::accumulate(deltaX.begin(), deltaX.begin() + elemIdx, 0.0)) * 2.0 / deltaX[elemIdx] - 1.0;

	ParamType map = 0.0;
	for (unsigned int i = 0; i < elemIdx; i++)
		map += deltaX[i];
	return  (x - map) * 2.0 / deltaX[elemIdx] - 1.0;
}
/**
 * @brief computes the Legendre-Gauss-Lobatto nodes and (inverse) quadrature weights
 * @param [in] polyDeg polynomial degree
 * @param [in, out] nodes Legendre Gauss Lobatto nodes
 * @param [in, out] invWeights Legendre Gauss quadrature weights
 * @param [in] invertWeights specifies if weights should be inverted
 */
void lglNodesWeights(const unsigned int polyDeg, Eigen::VectorXd& nodes, Eigen::VectorXd& invWeights, bool invertWeights = true);
/**
 * @brief computes the Legendre-Gauss nodes and quadrature weights
 * @detail Gauss quadrature exactly integrates polynomials up to degree 2N-1 with N points
 * @param [in] polyDeg polynomial degree (N-1 where N is number of points)
 * @param [in, out] nodes Legendre Gauss nodes
 * @param [in, out] weights Legendre Gauss quadrature weights
 * @param [in] invertWeights specifies if weights should be inverted
 */
void lgNodesWeights(const unsigned int polyDeg, Eigen::VectorXd& nodes, Eigen::VectorXd& weights, bool invertWeights = true);
/**
 * @brief computes the Legendre polynomial and its derivative
 * @param [in] polyDeg polynomial degree
 * @param [in, out] leg Legendre polynomial
 * @param [in, out] legDer Legendre polynomial derivative
 * @param [in] x evaluation point
 */
void legendrePolynomialAndDerivative(const int polyDeg, double& leg, double& legDer, const double x);
/**
 * @brief calculates the Gauss quadrature mass matrix for LGL interpolation polynomial on LG points
 * @detail exact integration of polynomials up to order 2N - 1
 * @param [in] LGLnodes Legendre Gauss Lobatto nodes
 * @param [in] nLGNodes number of Gauss quadrature nodes
 */
Eigen::MatrixXd gaussQuadratureMMatrix(const Eigen::VectorXd LGLnodes, const int nLGNodes);
/**
 * @brief calculates the barycentric weights for fast polynomial evaluation
 * @param [in] polyDeg polynomial degree
 * @param [in, out] baryWeights vector to store barycentric weights. Must already be initialized with ones!
 */
Eigen::VectorXd barycentricWeights(const unsigned int polyDeg, const Eigen::VectorXd nodes);
/**
 * @brief calculates the nodal (lagrange) polynomial derivative matrix
 * @param [in] polyDeg polynomial degree
 * @param [in] nodes polynomial interpolation nodes
 */
Eigen::MatrixXd derivativeMatrix(const unsigned int polyDeg, const Eigen::VectorXd nodes);
/**
 * @brief calculates the inverse mass matrix via transformation to orthonormal Jacobi (modal) basis
 * @detail the mass matrix used to compute integrals of the form \int_E \ell_i(\xi) \ell_j(\xi) (1 - \xi)^\alpha (1 + \xi)^\beta d\xi
 * @param [in] polyDeg polynomial degree
 * @param [in] nodes polynomial interpolation nodes
 * @param [in] alpha Jacobi polynomial coefficient
 * @param [in] beta Jacobi polynomial coefficient
 */
Eigen::MatrixXd invMMatrix(const unsigned int polyDeg, const Eigen::VectorXd nodes, const double alpha = 0.0, const double beta = 0.0);
/**
 * @brief calculates the mass matrix via transformation to orthonormal Jacobi (modal) basis
 * @detail mass matrix used to compute integrals of the form \int_E \ell_i(\xi) \ell_j(\xi) (1 - \xi)^\alpha (1 + \xi)^\beta d\xi
 * @param [in] polyDeg polynomial degree
 * @param [in] nodes polynomial interpolation nodes
 */
Eigen::MatrixXd mMatrix(const unsigned int polyDeg, const Eigen::VectorXd nodes, const double alpha, const double beta);
/**
 * @brief calculates the weighted mass matrix M^{(0,1)} for radial DG integrals
 * @detail For integrals of the form \int_E \ell_i(\xi) \ell_j(\xi) (1 + \xi) d\xi
 *         Used to construct radial weighted mass matrix: M_rho = (delta_rho/2) * M^{(0,1)} + rho_i * M^{(0,0)}
 * @param [in] polyDeg polynomial degree
 * @param [in] nodes polynomial interpolation nodes
 */
Eigen::MatrixXd weightedMMatrix(const unsigned int polyDeg, const Eigen::VectorXd nodes);
/**
 * @brief calculates a variable coefficient mass matrix via numerical quadrature
 * @detail Computes integrals of the form \int_E f(\xi) \ell_i(\xi) \ell_j(\xi) d\xi
 *         where f(\xi) is given by its values at LGL nodes.
 *         Uses Gauss-Legendre quadrature with nQuadPoints for exact integration.
 * @param [in] polyDeg polynomial degree
 * @param [in] LGLnodes LGL interpolation nodes
 * @param [in] coeffAtNodes coefficient values f evaluated at LGL nodes
 * @param [in] nQuadPoints number of Gauss quadrature points (default: polyDeg + 2 for safety)
 */
Eigen::MatrixXd varCoeffMMatrix(const unsigned int polyDeg, const Eigen::VectorXd LGLnodes,
                                 const Eigen::VectorXd coeffAtNodes, const int nQuadPoints = -1);
/**
 * @brief calculates a variable coefficient stiffness matrix via numerical quadrature
 * @detail Computes integrals of the form \int_E f(\xi) \ell'_i(\xi) \ell'_j(\xi) d\xi
 *         where f(\xi) is given by its values at LGL nodes.
 *         Used for dispersion terms with spatially varying D(x).
 * @param [in] polyDeg polynomial degree
 * @param [in] LGLnodes LGL interpolation nodes
 * @param [in] coeffAtNodes coefficient values f evaluated at LGL nodes
 * @param [in] nQuadPoints number of Gauss quadrature points (default: polyDeg + 2 for safety)
 */
Eigen::MatrixXd varCoeffStiffnessMatrix(const unsigned int polyDeg, const Eigen::VectorXd LGLnodes,
                                         const Eigen::VectorXd coeffAtNodes, const int nQuadPoints = -1);
/**
 * @brief calculates a mixed variable coefficient matrix for film diffusion terms
 * @detail Computes integrals of the form \int_E f(\xi) \ell_i(\xi) \ell_j(\xi) g(\xi) d\xi
 *         Used for radial film diffusion where k_f/rho terms appear.
 * @param [in] polyDeg polynomial degree
 * @param [in] LGLnodes LGL interpolation nodes
 * @param [in] coeff1AtNodes first coefficient values evaluated at LGL nodes
 * @param [in] coeff2AtNodes second coefficient values evaluated at LGL nodes (optional, default ones)
 * @param [in] nQuadPoints number of Gauss quadrature points (default: polyDeg + 3)
 */
Eigen::MatrixXd varCoeffMMatrixProduct(const unsigned int polyDeg, const Eigen::VectorXd LGLnodes,
                                        const Eigen::VectorXd coeff1AtNodes,
                                        const Eigen::VectorXd coeff2AtNodes = Eigen::VectorXd(),
                                        const int nQuadPoints = -1);
/**
 * @brief calculates a specific second order nodal stiffness matrix
 * @detail for integrals including terms of the form (1 - \xi)^\alpha (1 + \xi)^\beta. Computation via transformation to the respective Jacobi polynomial
 * @param [in] polyDeg polynomial degree
 * @param [in] a Jacobi polynomial parameter
 * @param [in] b Jacobi polynomial parameter
 * @param [in] nodes polynomial interpolation nodes
 */
Eigen::MatrixXd secondOrderStiffnessMatrix(const unsigned int polyDeg, const double alpha, const double beta, const Eigen::VectorXd nodes);
/**
 * @brief calculates the stiffness matrix via transformation to orthonormal Jacobi (modal) basis
 * @detail exact integration for integrals of the form \int_E \ell_i(\xi) \ell_j(\xi) (1 - \xi)^\alpha (1 + \xi)^\beta d\xi
 * @param [in] polyDeg polynomial degree
 * @param [in] nodes polynomial interpolation nodes
 * @param [in] alpha Jacobi polynomial coefficient
 * @param [in] beta Jacobi polynomial coefficient
 */
Eigen::MatrixXd stiffnessMatrix(const unsigned int polyDeg, const Eigen::VectorXd nodes, const double alpha = 0.0, const double beta = 0.0);
/**
 * @brief calculates the polynomial interpolation matrix between two sets of nodes
 * @param [in] newNodes set of nodes, the solution is interpolated to
 * @param [in] oldNodes set of nodes, the solution is interpolated from
 * @param [in] baryWeights barycentric weights of the polynomial to be interpolated
 */
Eigen::MatrixXd polynomialInterpolationMatrix(const Eigen::VectorXd newNodes, const Eigen::VectorXd oldNodes, const Eigen::VectorXd baryWeights);
/**
 * @brief returns a quadratic lifting matrix
 * @param [in] size quadratic matrix size
 */
Eigen::MatrixXd liftingMatrixQuadratic(const unsigned int size);
/**
 * @brief returns a (size x 2) lifting matrix
 * @param [in] rows number of matrix rows
 */
Eigen::MatrixXd liftingMatrix(const unsigned int size);
/**
 * @brief evaluates the jth Lagrange basis functions at given nodes
 * @param [in] j index of Lagrange basis function
 * @param [in] baseNodes interpolation nodes of Lagrange basis
 * @param [in] evalNodes evaluation nodes in [-1, 1]
 */
Eigen::VectorXd evalLagrangeBasis(const int j, const Eigen::VectorXd baseNodes, const Eigen::VectorXd evalNodes);
/**
 * @brief evaluates the derivative of the jth Lagrange basis function at given nodes
 * @param [in] j index of Lagrange basis function
 * @param [in] baseNodes interpolation nodes of Lagrange basis
 * @param [in] evalNodes evaluation nodes in [-1, 1]
 */
Eigen::VectorXd evalLagrangeBasisDerivative(const int j, const Eigen::VectorXd baseNodes, const Eigen::VectorXd evalNodes);
/**
 * @brief computes radial dispersion matrix S_g for a single cell with variable D(rho)
 * @detail For radial DG: S_g[i,j] = ∫ dL_i/dξ * L_j * ρ(ξ) * D(ρ(ξ)) dξ
 *         where ρ(ξ) = rho_left + (delta_rho/2) * (1 + ξ)
 * @param [in] polyDeg polynomial degree
 * @param [in] LGLnodes LGL interpolation nodes
 * @param [in] rho_left left boundary of cell in physical space
 * @param [in] delta_rho cell width in physical space
 * @param [in] dispAtNodes dispersion coefficient D evaluated at physical node positions
 * @param [in] nQuadPoints number of Gauss quadrature points (default: polyDeg + 2, use higher for nonlinear D)
 */
Eigen::MatrixXd radialDispersionMatrix(const unsigned int polyDeg, const Eigen::VectorXd LGLnodes,
                                        const double rho_left, const double delta_rho,
                                        const Eigen::VectorXd dispAtNodes, const int nQuadPoints = -1);
/**
 * @brief computes radial film diffusion mass matrix M_K for a single cell with variable k_f(rho)
 * @detail For radial DG: M_K[i,j] = ∫ L_i * L_j * ρ(ξ) * k_f(ρ(ξ)) dξ
 *         where ρ(ξ) = rho_left + (delta_rho/2) * (1 + ξ)
 * @param [in] polyDeg polynomial degree
 * @param [in] LGLnodes LGL interpolation nodes
 * @param [in] rho_left left boundary of cell in physical space
 * @param [in] delta_rho cell width in physical space
 * @param [in] kfAtNodes film diffusion coefficient k_f evaluated at physical node positions
 * @param [in] nQuadPoints number of Gauss quadrature points (default: polyDeg + 2, use higher for nonlinear k_f)
 */
Eigen::MatrixXd radialFilmDiffusionMatrix(const unsigned int polyDeg, const Eigen::VectorXd LGLnodes,
                                           const double rho_left, const double delta_rho,
                                           const Eigen::VectorXd kfAtNodes, const int nQuadPoints = -1);
/**
 * @brief evaluates the jth Lagrange basis functions at given nodes
 * @param [in, out] coords DG coordinate array
 * @param [in] nElem number of DG elements
 * @param [in] nNodes number of DG nodes (per element)
 * @param [in] DGnodes array of DG nodes
 * @param [in] length domain length
 * @param [in] leftElemBndries left element boundaries array
 */
void writeDGCoordinates(double* coords, const int nElem, const int nNodes, const double* DGnodes, const double length, const double* leftElemBndries);

} // namespace dgtoolbox
} // namespace parts
} // namespace model
} // namespace cadet

#endif  // LIBCADET_DGTOOLBOX_HPP_