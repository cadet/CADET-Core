(spline-interpolation)=

# Spline Interpolation

The spline interpolation model is a non-parametric, data-driven binding model.
It represents the equilibrium relation between pore-phase and solid-phase concentration
by interpolation of tabulated data.

Two interpolation modes are supported:

- **INDEPENDENT**: each component's equilibrium loading is a 1D cubic spline of its
  own pore-phase concentration.
- **COMPETITIVE_REGULAR_GRID**: all equilibrium loadings depend simultaneously on the full
  pore-phase concentration vector via multilinear interpolation on a regular
  Cartesian-product grid.

## INDEPENDENT mode

For each component $i$ and bound state $m$, an equilibrium loading

$$
c^{s,\ast}_{i,m} = f_{i,m}(c_{p,i})
$$

is constructed from user-provided data pairs $(\vec{c}^p_{i}, \vec{c}^s_{i,m})$.
Here, $c^p_{i}$ and $c^s_{i,m}$ denote the pore- and solid-phase concentration
of component $i$ and bound state $m$.

The spline function $f_{i,m}$ is a piecewise cubic polynomial. On an interval
$[c_{p,i}^{(k)}, c_{p,i}^{(k+1)}]$, it is evaluated as

$$
c^{s, \ast}_{i,m}(c_{p,i}) =
a_{i,m}^{(k)} h^3 + b_{i,m}^{(k)} h^2 + c_{i,m}^{(k)} h + d_{i,m}^{(k)},
\qquad
h = c_{p,i} - c_{p,i}^{(k)}.
$$

The spline coefficients are generated from the tabulated data using a cubic spline
construction. A monotonicity correction is applied to avoid non-physical oscillations
between supporting points.

### Extrapolation behavior

If $c^p_{i}$ lies outside the tabulated concentration range, the model extrapolates
using mixed boundary conditions:

- at $c^p_i < \operatorname{min}\left(\vec{c}^p_{i}\right)$, the second derivative is set to zero, resulting in a linear continuation beyond the boundary,
- at $c^p_i > \operatorname{max}\left(\vec{c}^p_{i}\right)$, the first derivative is set to zero, resulting in a flat continuation beyond the boundary.

## COMPETITIVE_REGULAR_GRID mode

In competitive mode, every bound state's equilibrium loading depends on the complete
pore-phase concentration vector
$\boldsymbol{c}_p = (c_{p,1}, \ldots, c_{p,N_c})$:

$$
c^{s,\ast}_{i,m} = f_{i,m}(\boldsymbol{c}_p).
$$

The user supplies a sorted support axis
$G_j = \{c_{p,j}^{(0)}, \ldots, c_{p,j}^{(K_j)}\}$ for each component $j$.
The evaluation grid is the Cartesian product $G_1 \times \cdots \times G_{N_c}$,
and equilibrium values $c^{s,\ast}_{i,m}$ must be provided at every one of the
$\prod_j (K_j + 1)$ grid points.

For a query point $\boldsymbol{c}_p$ located in the hypercell
$\prod_j [c_{p,j}^{(k_j)}, c_{p,j}^{(k_j+1)}]$, the interpolated value is the
multilinear (tensor-product linear) combination over the $2^{N_c}$ cell corners:

$$
f_{i,m}(\boldsymbol{c}_p)
=
\sum_{\boldsymbol{s} \in \{0,1\}^{N_c}}
\left(\prod_{j=1}^{N_c} w_j^{(s_j)}\right)
c^{s,\ast}_{i,m}\!\left(c_{p,j}^{(k_j+s_j)}\right)_{j=1}^{N_c},
$$

where the interpolation weights are
$w_j^{(0)} = 1 - \xi_j$, $w_j^{(1)} = \xi_j$, and

$$
\xi_j = \frac{c_{p,j} - c_{p,j}^{(k_j)}}{c_{p,j}^{(k_j+1)} - c_{p,j}^{(k_j)}}.
$$

This is equivalent to `scipy.interpolate.RegularGridInterpolator` with
`method='linear'`.

### Extrapolation behavior

If any pore-phase concentration lies outside the tabulated range, the multilinear
polynomial of the clamped boundary cell is extended linearly in each out-of-bounds
dimension. No constant fill value is applied.

## Porosity factor

In both modes, the interpolated value is scaled by an optional porosity factor
$\beta_{\text{poros},i,m}$ (`POROSITY_FACTOR`),

$$
c^{s,\ast}_{i,m} = \beta_{\text{poros},i,m} f_{i,m}(\boldsymbol{c}_p).
$$

Since equilibrium loadings are often reported on different reference volumes or masses,
this factor can be used as a simple correction for differences in porosity or unit
conversions between the tabulated data and the simulation conditions.
It may be given as a single value, which is shared by all bound states, or as one value
per bound state in component-major ordering, and defaults to $1$.

Note that, unlike the {ref}`neural-network` and {ref}`gaussian-process-regression`
binding models, no zero-concentration offset is subtracted here: the interpolant
reproduces the tabulated data exactly, so $c^{s,\ast}_{i,m}(0) = 0$ is obtained
simply by including the origin among the support points.

## Kinetic form

The model is used in a kinetic linear-driving-force form for both modes.
For each component $i$ and bound state $m$, the exchange term is based on
the deviation of the current solid-phase loading $c^s_{i,m}$ from the
interpolated equilibrium loading $c^{s,\ast}_{i,m}$:

$$
\frac{\partial c^s_{i,m}}{\partial t}
=
k^{\mathrm{kin}}_{i,m}\left(c^{s, \ast}_{i,m} - c^s_{i,m}\right).
$$

Thus, the spline interpolation model provides the equilibrium target, while the kinetic
constant $k^{\mathrm{kin}}_{i,m}$ controls how fast the equilibrium state is approached.

## Parameter sensitivities

Sensitivities for the spline interpolation data $\vec{c}^p_{i}, \vec{c}^s_{i,m}$ are currently not available.
Sensitivities for the kinetic constant $k^{\mathrm{kin}}_{i,m}$ are enabled.

## CADET-Core interface

For more information on model parameters required to define in CADET file format, see {ref}`spline-interpolation-config`.
