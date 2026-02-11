.. _hic_unified_config:

HIC Unified
~~~~~~~~~~~

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = HIC_UNIFIED**

For information on model equations, refer to :ref:`hic_unified_model`.


``IS_KINETIC``
   Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
   quasi-stationary. If a single value is given, the mode is set for all
   bound states. Otherwise, the adsorption mode is set for each bound
   state separately.

===================  =========================  =======================
**Type:** int        **Range:** {0,1}           **Length:** 1/NTOTALBND
===================  =========================  =======================

``HICUNI_KA``
   Adsorption rate constant

**Unit:** :math:`m_{MP}^{3}~m_{SP}^{-3}~s^{-1}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================

``HICUNI_KA_LIN``
   Coefficient of linear dependence of adsorption rate constant on
   modifier component

**Unit:** :math:`\text{[Mod]}^{-1}`

===================  =========================  
**Type:** double     **Length:** NCOMP
===================  ========================= 

``HICUNI_KD``
   Desorption rate constant

**Unit:** :math:`s^{-1}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================


``HICUNI_KS``
   Salt coefficient of adsorption; difference of
   water-protein and salt-protein interactions

===================  =========================
**Type:** double     **Length:** NCOMP
===================  =========================

``HICUNI_KP``
   Protein coefficient of adsorption; difference of
   water-protein and protein-protein interactions

**Unit:** :math:`m_{MP}^{3} mol^{-1}`

===================  =========================
**Type:** double     **Length:** NCOMP
===================  =========================

``HICUNI_NU``
   Number of ligands per ligand-protein interaction

**Unit:** [-]

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================


``HICUNI_NU_LIN``
   Coefficient of linear dependence of number of ligands per ligand-protein interaction
   on modifier component

**Unit:** :math:`\text{[Mod]}^{-1}`

===================  =========================
**Type:** double     **Length:** NCOMP
===================  =========================


``HICUNI_EPSILON``
   Empirical coefficient of interaction between bound molecules

**Unit:** :math:`m_{SP}^{3}~mol^{-1}`

===================  =========================
**Type:** double     **Length:** NCOMP
===================  =========================



``HICUNI_QMAX``
   Maximum binding capacity

**Unit:** :math:`mol~m_{SP}^{-3}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================


``HICUNI_BETA0``
   Parameters describing the number of highly ordered water molecules 
   that stabilize the hydrophobic surfaces at infinitely diluted 
   salt concentration

**Unit:** [-]

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** 1
===================  =========================  =========================================

``HICUNI_BETA1``
   Parameters describing the change in the number of highly ordered  
   water molecules that stabilize the hydrophobic surfaces with
   respect to changes in the salt concentration

**Unit:** :math:`m_{MP}^{3}~mol^{-1}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** 1
===================  =========================  =========================================

``HICUNI_REFC1``
   Reference for component 1 (optional, defaults to :math:`0.0`)

**Unit:** :math:`\text{[Mod]}`

===================   =========================================
**Type:** double      **Length:** 1
===================   =========================================

``HICUNI_RHO``
   Osmotic effect of :math:`c_0` on the water activity,
   calculated as osmotic_coefficient * molar_weight_of_water * ion_number
   (optional, defaults to :math:`3.35 * 10^{-5}`)

**Unit:** :math:`m_{MP}^{3}~mol^{-1}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\gt 0`   **Length:** 1
===================  =========================  =========================================
