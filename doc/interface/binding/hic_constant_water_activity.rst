.. _hic_constant_water_activity_config:

HIC Constant Water Activity
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = HIC_CONSTANT_WATER_ACTIVITY**

For information on model equations, refer to :ref:`hic_constant_water_activity_model`.


``IS_KINETIC``
   Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
   quasi-stationary. If a single value is given, the mode is set for all
   bound states. Otherwise, the adsorption mode is set for each bound
   state separately.

===================  =========================  =======================
**Type:** int        **Range:** {0,1}           **Length:** 1/NTOTALBND
===================  =========================  =======================

``HICCWA_KA``
   Adsorption rate constant

**Unit:** :math:`m_{MP}^{3}~m_{SP}^{-3}~s^{-1}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================

``HICCWA_KD``
   Desorption rate constant

**Unit:** :math:`s^{-1}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================

``HICCWA_NU``
   Number of ligands per ligand-protein interaction

**Unit:** [-]

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================


``HICCWA_QMAX``
   Maximum binding capacity

**Unit:** :math:`mol~m_{SP}^{-3}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================


``HICCWA_BETA0``
   Parameters describing the number of highly ordered water molecules 
   that stabilize the hydrophobic surfaces at infinitely diluted 
   salt concentration

**Unit:** [-]

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** 1
===================  =========================  =========================================

``HICCWA_BETA1``
   Parameters describing the change in the number of highly ordered  
   water molecules that stabilize the hydrophobic surfaces with
   respect to changes in the salt concentration

**Unit:** :math:`m_{MP}^{3}~mol^{-1}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** 1
===================  =========================  =========================================

