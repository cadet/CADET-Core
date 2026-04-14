.. _colloidal_particle_adsorption_config:

Colloidal Particle Adsorption
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = COLLOIDAL_PARTICLE_ADSORPTION**

For information on model equations, refer to :ref:`colloidal_particle_adsorption_model`.


``IS_KINETIC``
   So far, only "kinetic" has been implemented.
   Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
   quasi-stationary. If a single value is given, the mode is set for all
   bound states. Otherwise, the adsorption mode is set for each bound
   state separately.

===================  =========================  =========================================
**Type:** int        **Range:** {1}           **Length:** 1/NTOTALBND
===================  =========================  =========================================

``CPA_TEMPERATURE``
   Absolute temperature :math:`T`

**Unit:** :math:`\mathrm{K}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\gt 0`   **Length:** 1
===================  =========================  =========================================

``CPA_IONIC_STRENGTH``
   Ionic strength :math:`I_m` of the mobile phase. Ignored if
   ``CPA_SALT_IDX`` is set (ionic strength is then read from the
   corresponding mobile phase component).

**Unit:** :math:`\mathrm{mol \, m^{-3}}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** 1
===================  =========================  =========================================

``CPA_PERMITTIVITY``
   Relative permittivity :math:`\varepsilon` of the solvent

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\gt 0`   **Length:** 1
===================  =========================  =========================================

``CPA_SURFACE_DENSITY``
   Ligand surface density :math:`\Gamma_L`

**Unit:** :math:`\mathrm{mol \, m^{-2}}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** 1
===================  =========================  =========================================

``CPA_CHARGE_FULL_LIGAND``
   Charge of the fully protonated ligand :math:`\zeta_L`

===================  =========================  =========================================
**Type:** double                                **Length:** 1
===================  =========================  =========================================

``CPA_PK_LIGAND``
   Dissociation constant :math:`\mathrm{p}K_L` of the ligand

===================  =========================  =========================================
**Type:** double                                **Length:** 1
===================  =========================  =========================================

``CPA_PH_REF``
   Reference pH value :math:`\mathrm{pH}_{\mathrm{ref}}` for the
   quadratic protein charge model

===================  =========================  =========================================
**Type:** double                                **Length:** 1
===================  =========================  =========================================

``CPA_SURFACE_AREA``
   Specific adsorber surface area per skeleton volume :math:`A_{s,i}`

**Unit:** :math:`\mathrm{m^{-1}}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\gt 0`   **Length:** NCOMP
===================  =========================  =========================================

``CPA_PROTEIN_RADIUS``
   Hydrodynamic radius :math:`a_i` of each protein component

**Unit:** :math:`\mathrm{m}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================

``CPA_COMP_CHARGE_REF``
   Reference protein net charge :math:`Z_{i,\mathrm{ref}}` at
   :math:`\mathrm{pH}_{\mathrm{ref}}`

===================  =========================  =========================================
**Type:** double                                **Length:** NCOMP
===================  =========================  =========================================

``CPA_COMP_CHARGE_LIN``
   Linear coefficient :math:`Z_{i,\mathrm{lin}}` of the
   pH-dependent protein net charge

===================  =========================  =========================================
**Type:** double                                **Length:** NCOMP
===================  =========================  =========================================

``CPA_COMP_CHARGE_QUAD``
   Quadratic coefficient :math:`Z_{i,\mathrm{quad}}` of the
   pH-dependent protein net charge

===================  =========================  =========================================
**Type:** double                                **Length:** NCOMP
===================  =========================  =========================================

``CPA_COMP_LAT_CHARGE``
   Lateral charge :math:`Z_{\mathrm{lat},i}` used for computing
   the pairwise Yukawa coefficient :math:`\beta_{ij}` in the lateral
   protein–protein interaction

===================  =========================  =========================================
**Type:** double                                **Length:** NCOMP
===================  =========================  =========================================

``CPA_DELTA_REF``
   Reference value :math:`\delta_{i,\mathrm{ref}}` for the
   logarithmic interaction layer thickness parameterisation

===================  =========================  =========================================
**Type:** double                                **Length:** NCOMP
===================  =========================  =========================================

``CPA_DELTA_LIN``
   Linear coefficient :math:`\delta_{i,\mathrm{lin}}` for the
   logarithmic interaction layer thickness parameterisation

===================  =========================  =========================================
**Type:** double                                **Length:** NCOMP
===================  =========================  =========================================

``CPA_DIFFUSION_COEFF``
   Pore diffusion coefficient :math:`D_i` used in the kinetic rate
   constant :math:`k_{\mathrm{kin},i}`

**Unit:** :math:`\mathrm{m^{2} \, s^{-1}}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================

``CPA_PROTON_IDX``
   0-based index of the proton/pH component in the mobile phase
   (optional, defaults to 0). The corresponding component must be
   non-binding (``NBOUND = 0``).

===================  =========================  =========================================
**Type:** int        **Range:** :math:`\ge 0`   **Length:** 1
===================  =========================  =========================================

``CPA_SALT_IDX``
   0-based index of the salt component in the mobile phase (optional).
   If set, the ionic strength :math:`I_m` is read from this mobile phase
   component instead of from ``CPA_IONIC_STRENGTH``. The corresponding
   component must be non-binding (``NBOUND = 0``).

===================  =========================  =========================================
**Type:** int        **Range:** :math:`\ge 0`   **Length:** 1
===================  =========================  =========================================

``CPA_MAXITER``
   Maximum number of Newton iterations for solving the adsorber surface
   potential :math:`\psi_{0,A}` (optional, defaults to 100)

===================  =========================  =========================================
**Type:** int        **Range:** :math:`\ge 1`   **Length:** 1
===================  =========================  =========================================
