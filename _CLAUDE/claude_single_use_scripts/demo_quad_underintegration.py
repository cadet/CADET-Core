"""Demonstrate the under-integration caused by the old quadrature node count.

Builds the dispersion- and area-weighted mass matrix

    M[i,j] = int_{-1}^{1} l_i(xi) l_j(xi) A(xi) D(xi) dxi

for a representative frustum case (polyDeg = 3, DISPERSION_SPATIAL_DEPENDENCE_POLYDEG = 2),
whose integrand has degree 2*polyDeg + geomFactorDegree + axDispQuadDeg = 10, and compares
the old node count (5, exact only up to degree 9) and the fixed one (6) against a
converged reference.
"""

import numpy as np

POLYDEG = 3
AX_DISP_QUAD_DEG = 2
GEOM_DEGREE = 2  # frustum, A ~ r(x)^2

degree = AX_DISP_QUAD_DEG + GEOM_DEGREE + 2 * POLYDEG
old_nodes = (AX_DISP_QUAD_DEG + GEOM_DEGREE + 2 * POLYDEG + 1) // 2
new_nodes = (AX_DISP_QUAD_DEG + GEOM_DEGREE + 2 * POLYDEG + 2) // 2
print(f"integrand degree = {degree}, old node count = {old_nodes} "
      f"(exact to degree {2 * old_nodes - 1}), fixed node count = {new_nodes} "
      f"(exact to degree {2 * new_nodes - 1})")


def lgl_nodes(poly_deg):
    """Legendre-Gauss-Lobatto nodes on [-1, 1]."""
    if poly_deg == 1:
        return np.array([-1.0, 1.0])
    # interior nodes are the roots of P'_n
    coeffs = np.zeros(poly_deg + 1)
    coeffs[poly_deg] = 1.0
    interior = np.polynomial.legendre.Legendre(coeffs).deriv().roots()
    return np.concatenate(([-1.0], interior, [1.0]))


def lagrange(j, base, x):
    out = np.ones_like(np.asarray(x, dtype=float))
    for m in range(len(base)):
        if m != j:
            out = out * (x - base[m]) / (base[j] - base[m])
    return out


base = lgl_nodes(POLYDEG)

# area (degree 2) and dispersion dependence (degree 2) with arbitrary but generic coefficients
area = np.polynomial.Polynomial([1.0, 0.4, 0.15])
disp = np.polynomial.Polynomial([1.0, -0.3, 0.2])


def mass_matrix(n_quad):
    nodes, weights = np.polynomial.legendre.leggauss(n_quad)
    basis = np.array([lagrange(i, base, nodes) for i in range(len(base))])
    weight = weights * area(nodes) * disp(nodes)
    return (basis * weight) @ basis.T


reference = mass_matrix(40)
for n in (old_nodes, new_nodes):
    err = np.max(np.abs(mass_matrix(n) - reference)) / np.max(np.abs(reference))
    label = "old" if n == old_nodes else "fixed"
    print(f"{label:5s} ({n} nodes): max relative deviation from the exact integral = {err:.3e}")
