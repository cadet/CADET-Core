.. _spline_interpolation:

Spline Interpolation
~~~~~~~~~~~~~~~~~~~~~

The spline interpolation model is a non-parametric, data-driven binding model.
It represents the equilibrium relation between pore-phase and solid-phase concentration
by interpolation of tabulated data.

Two interpolation modes are supported:

* **INDEPENDENT**: each component's equilibrium loading is a 1D cubic spline of its
  own pore-phase concentration.
* **COMPETITIVE_REGULAR_GRID**: all equilibrium loadings depend simultaneously on the full
  pore-phase concentration vector via multilinear interpolation on a regular
  Cartesian-product grid.

INDEPENDENT mode
*****************

For each component :math:`i` and bound state :math:`m`, an equilibrium loading

.. math::

    c^{s,\ast}_{i,m} = f_{i,m}(c_{p,i})

is constructed from user-provided data pairs :math:`(\vec{c}^p_{i}, \vec{c}^s_{i,m})`.
Here, :math:`c^p_{i}` and :math:`c^s_{i,m}` denote the pore- and solid-phase concentration
of component :math:`i` and bound state :math:`m`.

The spline function :math:`f_{i,m}` is a piecewise cubic polynomial. On an interval
:math:`[c_{p,i}^{(k)}, c_{p,i}^{(k+1)}]`, it is evaluated as

.. math::

    c^{s, \ast}_{i,m}(c_{p,i}) =
    a_{i,m}^{(k)} h^3 + b_{i,m}^{(k)} h^2 + c_{i,m}^{(k)} h + d_{i,m}^{(k)},
    \qquad
    h = c_{p,i} - c_{p,i}^{(k)}.

The spline coefficients are generated from the tabulated data using a cubic spline
construction. A monotonicity correction is applied to avoid non-physical oscillations
between supporting points.

Extrapolation behavior
-----------------------

If :math:`c^p_{i}` lies outside the tabulated concentration range, the model extrapolates
using mixed boundary conditions:

- at :math:`c^p_i < \operatorname{min}\left(\vec{c}^p_{i}\right)`, the second derivative is set to zero, resulting in a linear continuation beyond the boundary,
- at :math:`c^p_i > \operatorname{max}\left(\vec{c}^p_{i}\right)`, the first derivative is set to zero, resulting in a flat continuation beyond the boundary.

COMPETITIVE_REGULAR_GRID mode
*****************************

In competitive mode, every bound state's equilibrium loading depends on the complete
pore-phase concentration vector
:math:`\boldsymbol{c}_p = (c_{p,1}, \ldots, c_{p,N_c})`:

.. math::

    c^{s,\ast}_{i,m} = f_{i,m}(\boldsymbol{c}_p).

The user supplies a sorted support axis
:math:`G_j = \{c_{p,j}^{(0)}, \ldots, c_{p,j}^{(K_j)}\}` for each component :math:`j`.
The evaluation grid is the Cartesian product :math:`G_1 \times \cdots \times G_{N_c}`,
and equilibrium values :math:`c^{s,\ast}_{i,m}` must be provided at every one of the
:math:`\prod_j (K_j + 1)` grid points.

For a query point :math:`\boldsymbol{c}_p` located in the hypercell
:math:`\prod_j [c_{p,j}^{(k_j)}, c_{p,j}^{(k_j+1)}]`, the interpolated value is the
multilinear (tensor-product linear) combination over the :math:`2^{N_c}` cell corners:

.. math::

    f_{i,m}(\boldsymbol{c}_p)
    =
    \sum_{\boldsymbol{s} \in \{0,1\}^{N_c}}
    \left(\prod_{j=1}^{N_c} w_j^{(s_j)}\right)
    c^{s,\ast}_{i,m}\!\left(c_{p,j}^{(k_j+s_j)}\right)_{j=1}^{N_c},

where the interpolation weights are
:math:`w_j^{(0)} = 1 - \xi_j`, :math:`w_j^{(1)} = \xi_j`, and

.. math::

    \xi_j = \frac{c_{p,j} - c_{p,j}^{(k_j)}}{c_{p,j}^{(k_j+1)} - c_{p,j}^{(k_j)}}.

This is equivalent to ``scipy.interpolate.RegularGridInterpolator`` with
``method='linear'``.

Extrapolation behavior
-----------------------

If any pore-phase concentration lies outside the tabulated range, the multilinear
polynomial of the clamped boundary cell is extended linearly in each out-of-bounds
dimension. No constant fill value is applied.

Kinetic form
*************

The model is used in a kinetic linear-driving-force form for both modes.
For each component :math:`i` and bound state :math:`m`, the exchange term is based on
the deviation of the current solid-phase loading :math:`c^s_{i,m}` from the
interpolated equilibrium loading :math:`c^{s,\ast}_{i,m}`:

.. math::

    \frac{\partial c^s_{i,m}}{\partial t}
    =
    k^{\mathrm{kin}}_{i,m}\left(c^{s, \ast}_{i,m} - c^s_{i,m}\right).

Thus, the spline interpolation model provides the equilibrium target, while the kinetic
constant :math:`k^{\mathrm{kin}}_{i,m}` controls how fast the equilibrium state is approached.

Parameter sensitivities
***********************

Sensitivities for the spline interpolation data :math:`\vec{c}^p_{i}, \vec{c}^s_{i,m}` are currently not available.
Sensitivities for the kinetic constant :math:`k^{\mathrm{kin}}_{i,m}` are enabled.

CADET-Core interface
********************

For more information on model parameters required to define in CADET file format, see :ref:`spline_interpolation_config`.