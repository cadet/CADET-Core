.. _multi_component_colloidal_config:

Multi Component Colloidal
~~~~~~~~~~~~~~~~~~~~~~~~~

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = MULTI_COMPONENT_COLLOIDAL**


``IS_KINETIC``
   Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
   quasi-stationary. If a single value is given, the mode is set for all
   bound states. Otherwise, the adsorption mode is set for each bound
   state separately.

===================  =========================  =========================================
**Type:** int        **Range:** {0,1}  		    **Length:** 1/NTOTALBND
===================  =========================  =========================================

``COL_PHI``
   Phase ratio :math:`\Phi`

**Unit:** :math:`m^{2} m_{s}^{-3}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** 1
===================  =========================  =========================================

``COL_KAPPA_EXP``
   Screening term exponent factor :math:`\kappa_{ef}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** 1
===================  =========================  =========================================

``COL_KAPPA_FACT``
   Screening term factor :math:`\kappa_{f}`

**Unit:** :math:`m \cdot mM^{-1}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** 1
===================  =========================  =========================================

``COL_KAPPA_CONST``
   Screening term constant :math:`\kappa_{c}`

**Unit:** :math:`m`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** 1
===================  =========================  =========================================

``COL_CORDNUM``
   Coordination number :math:`n`

===================  =========================  =========================================
**Type:** int        **Range:** :math:`\ge 0`   **Length:** 1
===================  =========================  =========================================

``COL_LOGKEQ_PH_EXP``
   Protein-resin interaction :math:`K_{e,i}` equilibrium: Constant exponent :math:`k_{e,i}` for pH
   If pH is not considered, this value will be not be used but must still be specified, i.e. can be any number.

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NTOTALBND
===================  =========================  =========================================

``COL_LOGKEQ_SALT_POWEXP``
   Protein-resin interaction :math:`K_{e,i}` equilibrium: Constant pre-factor :math:`k_{a,i}` for salt power

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NTOTALBND
===================  =========================  =========================================

``COL_LOGKEQ_SALT_POWFAC``
   Protein-resin interaction :math:`K_{e,i}` equilibrium: Constant exponent :math:`k_{b,i}` for salt power

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NTOTALBND
===================  =========================  =========================================

``COL_LOGKEQ_SALT_EXPFAC``
   Protein-resin interaction :math:`K_{e,i}` equilibrium: Constant pre-factor :math:`k_{c,i}` for e-function with salt power

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NTOTALBND
===================  =========================  =========================================

``COL_LOGKEQ_SALT_EXPARGMULT``
   Protein-resin interaction :math:`K_{e,i}` equilibrium: Constant power factor :math:`k_{d,i}` for salt in e-function

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NTOTALBND
===================  =========================  =========================================

``COL_BPP_PH_EXP``
   Protein-protein interaction :math:`b_{pp,i}`: Constant power term :math:`b_{e,i}` for pH.
   If pH is not considered, this value will be not be used but must still be specified, i.e. can be any number.

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NTOTALBND
===================  =========================  =========================================

``COL_BPP_SALT_POWFACT``
   Protein-protein interaction :math:`b_{pp,i}`: Constant power pre-factor :math:`b_{a,i}` for salt

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NTOTALBND
===================  =========================  =========================================

``COL_BPP_SALT_POWEX``
   Protein-protein interaction :math:`b_{pp,i}`: Constant power :math:`b_{b,i}` for salt

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NTOTALBND
===================  =========================  =========================================

``COL_BPP_SALT_EXPFACT``
   Protein-protein interaction :math:`b_{pp,i}`: Constant pre-factor :math:`b_{c,i}` e-function with salt power

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NTOTALBND
===================  =========================  =========================================

``COL_BPP_SALT_EXPARGMULT``
   Protein-protein interaction :math:`b_{pp,i}`: Constant power factor :math:`b_{d,i}` for salt in e-function

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NTOTALBND
===================  =========================  =========================================

``COL_PROTEIN_RADIUS``
   Protein radius :math:`r_i`

**Unit:** :math:`m`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NTOTALBND
===================  =========================  =========================================

``COL_KKIN``
   Adsorption rate constants :math:`K_\text{kin}` in state-major ordering

**Unit:** :math:`s^{-1}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NTOTALBND
===================  =========================  =========================================

``COL_LINEAR_THRESHOLD``
   Threshold concentration :math:`c_\epsilon` for switching to linear approximation

===================  =========================  =========================================
**Type:** double     **Range:** :math:`> 0`     **Length:** 1
===================  =========================  =========================================

``COL_USE_PH``
   Selects if pH is included in the model or not: 1 = yes, 0 = no.

===================  =========================  =========================================
**Type:** int        **Range:** {0,1}           **Length:** 1
===================  =========================  =========================================
