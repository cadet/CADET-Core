.. _affinity_complex_titration_config:

Affinity Complex Titration
==========================

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = AFFINITY_COMPLEX_TITRATION**

For information on model equations, refer to :ref:`affinity_complex_titration`. The first component should be pH. 


``IS_KINETIC``
   Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
   quasi-stationary. If a single value is given, the mode is set for all
   bound states. Otherwise, the adsorption mode is set for each bound
   state separately.

===================  =========================  =========================================
**Type:** int        **Range:** {0,1}  		    **Length:** 1/NTOTALBND
===================  =========================  =========================================

``ACT_KA``
   Adsorption rate constants

**Unit:** :math:`m_{MP}^3~mol^{-1}~s^{-1}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================

``ACT_KD``
   Desorption rate constants

**Unit:** :math:`s^{-1}`

===================  =========================  ==================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  ================================== 

``ACT_QMAX``
   Maximum binding capacities

**Unit:** :math:`mol~m_{SP}^{-3}`

===================  =========================  ==================================
**Type:** double     **Range:** :math:`\gt 0`   **Length:** NCOMP
===================  =========================  ================================== 

``ACT_ETAA``
   Hill-type coefficients denoting the slope for the binding capacity changes as a function of pH changes

**Unit:** :math:`1`

===================  =============================  ==================================
**Type:** double     **Range:** :math:`\mathbb{R}`   **Length:** NCOMP
===================  =============================  ================================== 

``ACT_ETAG``
   Hill-type coefficients denoting the slope for the equilibrium constant changes as a function of pH changes

**Unit:** :math:`1`

===================  =============================  ==================================
**Type:** double     **Range:** :math:`\mathbb{R}`   **Length:** NCOMP
===================  =============================  ================================== 

``ACT_PKAA``
   Center point for the binding capacity changes as a function of pH changes

**Unit:** :math:`1`

===================  ==============================  ==================================
**Type:** double     **Range:** :math:`\mathbb{R}`   **Length:** NCOMP
===================  ==============================  ================================== 

``ACT_PKAG``
   Center point for the equilibrium constant changes as a function of pH changes

**Unit:** :math:`1`

===================  ==============================  ==================================
**Type:** double     **Range:** :math:`\mathbb{R}`   **Length:** NCOMP
===================  ==============================  ================================== 