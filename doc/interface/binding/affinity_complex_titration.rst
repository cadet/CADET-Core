.. _affinity_complex_titration_config:

Affinity Complex Titration
==========================

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = AFFINITY_COMPLEX_TITRATION**

For information on model equations, refer to :ref:`affinity_complex_titration`. 
The first component is either ion concentration (:math:`c_{\mathrm{H}^+}`, :math:`c_{\mathrm{Na}^+}`) or negative log ion concentration (:math:`\mathrm{pH}`, :math:`\mathrm{pNa}`).

``IS_KINETIC``
   Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
   quasi-stationary. If a single value is given, the mode is set for all
   bound states. Otherwise, the adsorption mode is set for each bound
   state separately.

===================  =========================  =========================================
**Type:** int        **Range:** {0,1}  		    **Length:** 1/NTOTALBND
===================  =========================  =========================================

``ACT_USE_SALT_CONC``
   Optional. Selects to use either ion concentration or negative log ion concentration as the first component. 

   * ``False``: input is negative log ion concentration :math:`\mathrm{pIon}`; use ``ACT_PKAA`` / ``ACT_PKAG``.
   * ``True``: input is ion concentration :math:`c_{\mathrm{Ion}^{n+}}`; use ``ACT_CMID_A`` / ``ACT_CMID_G``.

   Default is False. 

===================  =========================  =========================================
**Type:** bool       **Range:** {False, True}   **Length:** 1
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
   Hill-type coefficients denoting the slope for the binding capacity changes as a function of pH/pIon changes

**Unit:** :math:`1`

===================  =============================  ==================================
**Type:** double     **Range:** :math:`\mathbb{R}`   **Length:** NCOMP
===================  =============================  ================================== 

``ACT_ETAG``
   Hill-type coefficients denoting the slope for the equilibrium constant changes as a function of pH/pIon changes

**Unit:** :math:`1`

===================  =============================  ==================================
**Type:** double     **Range:** :math:`\mathbb{R}`   **Length:** NCOMP
===================  =============================  ================================== 

``ACT_PKAA``
   Center point for the binding capacity changes as a function of pH/pIon changes

**Unit:** :math:`1`

===================  ==============================  ==================================
**Type:** double     **Range:** :math:`\mathbb{R}`   **Length:** NCOMP
===================  ==============================  ================================== 

``ACT_PKAG``
   Center point for the equilibrium constant changes as a function of pH/pIon changes

**Unit:** :math:`1`

===================  ==============================  ==================================
**Type:** double     **Range:** :math:`\mathbb{R}`   **Length:** NCOMP
===================  ==============================  ================================== 

``ACT_CMID_A``
   Midpoint ion concentration for the binding capacity changes. 
   Only used when ``ACT_USE_ION_CONC = True``.

**Unit:** :math:`mol~m_{MP}^{-3}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================

``ACT_CMID_G``
   Midpoint ion concentration for the equilibrium constant changes.
   Only used when ``ACT_USE_ION_CONC = True``.

**Unit:** :math:`mol~m_{MP}^{-3}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================