.. _kumar_langmuir_config:

Kumar-Langmuir
~~~~~~~~~~~~~~

**Group /input/model/unit_XXX/adsorption – ADSORPTION_MODEL = KUMAR_MULTI_COMPONENT_LANGMUIR**


``IS_KINETIC``
   Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
   quasi-stationary. If a single value is given, the mode is set for all
   bound states. Otherwise, the adsorption mode is set for each bound
   state separately.

===================  =========================  =========================================
**Type:** int        **Range:** {0,1}  		**Length:** 1/NTOTALBND
===================  =========================  =========================================

``KMCL_KA``
   Adsorption pre-exponential factors

**Unit:** :math:`m_{MP}^3~mol^{-1}~s^{-1}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================

``KMCL_KD``
   Desorption rate

**Unit:** :math:`m_{MP}^{3\nu_{i}}~mol^{-\nu_{i}}~s^{-1}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================

``KMCL_KACT``
   Activation temperatures

**Unit:** :math:`K`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================

``KMCL_QMAX``
   Maximum adsorption capacities

**Unit:** :math:`mol~m_{SP}^{-3}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\gt 0`   **Length:** NCOMP
===================  =========================  =========================================

``KMCL_NU``
   Salt exponents / characteristic charges

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\gt 0`   **Length:** NCOMP
===================  =========================  =========================================

``KMCL_TEMP``
   Temperature

**Unit:** :math:`K`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** 1
===================  =========================  =========================================
