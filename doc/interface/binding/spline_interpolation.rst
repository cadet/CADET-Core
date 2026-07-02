.. _spline_interpolation_config:

Spline Interpolation
~~~~~~~~~~~~~~~~~~~~~


**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = SPLINE_INTERPOLATION**

For information on model equations, refer to :ref:`spline_interpolation`.

``IS_KINETIC``
   Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
   quasi-stationary. If a single value is given, the mode is set for all
   bound states. Otherwise, the adsorption mode is set for each bound
   state separately.

===================  =========================  =======================
**Type:** int        **Range:** {0,1}           **Length:** 1/NTOTALBND
===================  =========================  =======================

``INTERPOLATION_MODE``
   Selects the interpolation mode:

   - ``INDEPENDENT``: each component's equilibrium loading is an independent
     1D cubic spline of its own pore-phase concentration.
   - ``COMPETITIVE``: all equilibrium loadings depend on the full pore-phase
     concentration vector via multilinear interpolation on a regular
     Cartesian-product grid (see :ref:`spline_interpolation`).

===================  ===================================
**Type:** string     **Value:** INDEPENDENT, COMPETITIVE
===================  ===================================

``CP_VALS_COMP_XXX``
   Pore-phase concentration support points for component ``XXX``.

   *INDEPENDENT mode*: the 1D spline knots in strictly ascending order.

   *COMPETITIVE mode*: sample values whose distinct entries define the axis
   for component ``XXX`` of the Cartesian-product grid. All components must
   supply the same number of samples, together covering the full grid (i.e.
   every combination of distinct per-component values must appear exactly once).

**Unit:** :math:`mol~m_{MP}^{-3}`

===================  =========================  =======================
**Type:** double     **Range:** :math:`\geq 0`  **Length:** N
===================  =========================  =======================

``CS_VALS_COMP_XXX_BND_YYY``
   Solid-phase equilibrium loadings corresponding to ``CP_VALS_COMP_XXX`` for bound state ``YYY`` of component ``XXX``.

**Unit:** :math:`mol~m_{SP}^{-3}`

===================  =========================  =======================
**Type:** double     **Range:** :math:`\geq 0`  **Length:** N
===================  =========================  =======================

``SPLINE_KKIN``
   Linear-driving-force coefficients in component-major ordering.

**Unit:** :math:`s^{-1}`

===================  =========================  =======================
**Type:** double     **Range:** :math:`\geq 0`  **Length:** NTOTALBND
===================  =========================  =======================
