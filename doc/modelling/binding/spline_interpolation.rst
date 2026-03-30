.. _spline_interpolation:

Spline Interpolation
~~~~~~~~~~~~~~~~~~~~~

The spline interpolation model is a non-parametric, data-driven binding model that represents the equilibrium relation between pore-phase and solid-phase concentration by cubic spline interpolation of tabulated data.

For each component :math:`i` and bound state :math:`m`, an equilibrium loading

.. math::

    c^{s,\ast}_{i,m} = f_{i,m}(c_{p,i})

is constructed from user-provided data pairs :math:`(\vec{c}^p_{i}, \vec{c}^s_{i,m})`.  
Here, :math:`c^p_{i},c^s_{i,m}` denote the pore- and solid-phase concentration of component :math:`i` and bound state :math:`m`, and :math:`c^{s,\ast}_{i,m}` is the corresponding equilibrium solid-phase loading.

The spline function :math:`f_{i,m}` is a piecewise cubic polynomial. On an interval
:math:`[c_{p,i}^{(k)}, c_{p,i}^{(k+1)}]`, it is evaluated as

.. math::

    c^{s, \ast}_{i,m}(c_{p,i}) =
    a_{i,m}^{(k)} h^3 + b_{i,m}^{(k)} h^2 + c_{i,m}^{(k)} h + d_{i,m}^{(k)},
    \qquad
    h = c_{p,i} - c_{p,i}^{(k)}.

The spline coefficients are generated from the tabulated data using a cubic spline construction.  
A monotonicity correction is applied to avoid non-physical oscillations between supporting points.

Kinetic form
************

The model is used in a kinetic linear-driving-force form.
For each component :math:`i` and bound state :math:`m`, the exchange term is based on the deviation of the current solid-phase loading :math:`c^s_{i,m}` from the spline-predicted equilibrium loading :math:`c^{s,\ast}_{i,m}`:

.. math::

    \frac{\partial c^s_{i,m}}{\partial t}
    =
    k^{\mathrm{kin}}_{i,m}\left(c^{s, \ast}_{i,m}(c^p_{i}) - c^s_{i,m}\right).

Thus, the spline interpolation model provides the equilibrium target, while the kinetic constant :math:`k^{\mathrm{kin}}_{i,m}` controls how fast the equilibrium state is approached.

Extrapolation behavior
**********************

If :math:`c^p_{i}` lies outside the tabulated concentration range, the model extrapolates using mixed boundary conditions:

- at :math:`c^p_i < \operatorname{min}\left(\vec{c}^p_{i}\right)`, the second derivative is set to zero, resulting in a linear continuation beyond the boundary,

- at :math:`c^p_i > \operatorname{max}\left(\vec{c}^p_{i}\right)`, the first derivative is set to zero, resulting in a flat continuation beyond the boundary.

Parameter sensitivities
***********************

Sensitivities for the spline interpolation data :math:`\vec{c}^p_{i}, \vec{c}^s_{i,m}` are currently not available.
Sensitivities for the kinetic constant :math:`k^{\mathrm{kin}}_{i,m}` are enabled.

CADET-Core interface
********************

For more information on model parameters required to define in CADET file format, see :ref:`spline_interpolation_config`.