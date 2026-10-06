.. _colloidal_particle_adsorption_config:

Colloidal Particle Adsorption
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = COLLOIDAL_PARTICLE_ADSORPTION**

For information on model equations, refer to :ref:`colloidal_particle_adsorption_model`.


``IS_KINETIC``
   Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
   quasi-stationary. If a single value is given, the mode is set for all
   bound states. Otherwise, the adsorption mode is set for each bound
   state separately.

===================  =========================  =========================================
**Type:** int        **Range:** {0,1}           **Length:** 1/NTOTALBND
===================  =========================  =========================================

``CPA_TEMPERATURE``
   Absolute temperature :math:`T`

**Unit:** :math:`\mathrm{K}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`> 0`     **Length:** 1
===================  =========================  =========================================

``CPA_IONIC_STRENGTH``
   Ionic strength :math:`I_m` of the mobile phase. Ignored if
   ``CPA_IONIC_VALENCE`` is set (ionic strength is then computed
   from the pore-phase concentrations).

**Unit:** :math:`\mathrm{mol \, m^{-3}}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`> 0`     **Length:** 1
===================  =========================  =========================================

``CPA_PERMITTIVITY``
   Relative permittivity :math:`\varepsilon` of the solvent

===================  =========================  =========================================
**Type:** double     **Range:** :math:`> 0`     **Length:** 1
===================  =========================  =========================================

``CPA_LIGAND_DENSITY``
   Ligand surface density :math:`\Gamma_L`

**Unit:** :math:`\mathrm{mol \, m^{-2}}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** 1
===================  =========================  =========================================

``CPA_LIGAND_CHARGE_FULL``
   Charge of the fully protonated ligand :math:`\zeta_L`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`     **Length:** 1
===================  =========================  =========================================

``CPA_LIGAND_PK``
   Dissociation constant :math:`\mathrm{p}K_L` of the ligand

===================  =========================  =========================================
**Type:** double     **Range:** :math:`> 0`     **Length:** 1
===================  =========================  =========================================

``CPA_PH_REF``
   Reference pH value :math:`\mathrm{pH}_{\mathrm{ref}}` for the
   protein charge polynomial

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`     **Length:** 1
===================  =========================  =========================================

``CPA_SPECIFIC_SURFACE_AREA``
   Specific adsorber surface area per skeleton volume :math:`A_{s,i}`.

**Unit:** :math:`\mathrm{m^{-1}}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`> 0`     **Length:** NCOMP
===================  =========================  =========================================

``CPA_RADIUS``
   Hydrodynamic radius :math:`a_i` of each protein component

**Unit:** :math:`\mathrm{m}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`> 0`     **Length:** NCOMP
===================  =========================  =========================================

``CPA_EFFECTIVE_CHARGE_COEF``
   Matrix of coefficients :math:`z_{i,k}` for the pH-dependent protein
   net charge. Matrix rows correspond to increasing polynomial powers
   :math:`k=0,\ldots,P`, and columns correspond to components. The first row
   therefore contains the constant charges :math:`z_{i,0}` at
   :math:`\mathrm{pH}=\mathrm{pH}_{\mathrm{ref}}`. The net charge is

   .. math::

      Z_i(\mathrm{pH}) = z_{i,0} + \sum_{k=1}^{P} z_{i,k}
      \left(\mathrm{pH}-\mathrm{pH}_{\mathrm{ref}}\right)^k.

   The matrix is supplied in polynomial-order-row-major ordering.
   Every row must contain one value for every component, including non-binding
   components.
   Coefficients that do not apply should be set to zero.

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`     **Length:** NCOMP * (P + 1)
===================  =========================  =========================================

``CPA_LAT_CHARGE``
   Lateral charge :math:`Z_{\mathrm{lat},i}` used for computing
   the pairwise Yukawa coefficient :math:`\beta_{ij}` in the lateral
   protein–protein interaction

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`     **Length:** NCOMP
===================  =========================  =========================================

``CPA_DELTA_REF``
   Positive reference interaction layer thickness :math:`\delta_{i,\mathrm{ref}}`.

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`     **Length:** NCOMP
===================  =========================  =========================================

``CPA_DELTA_LIN``
   Linear coefficient :math:`\delta_{i,\mathrm{lin}}` for the
   logarithmic interaction layer thickness parameterisation

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`      **Length:** NCOMP
===================  =========================  =========================================

``CPA_KKIN``
   Kinetic prefactor :math:`k^*_{\mathrm{kin},i}` used to calculate the
   kinetic rate constant :math:`k_{\mathrm{kin},i}`.

**Unit:** :math:`\mathrm{s^{-1}}`

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

``CPA_IONIC_VALENCE``
   Integer valence charge :math:`z_i` for each component (optional).
   If provided, the ionic strength is computed from the pore-phase
   concentrations, and the Davies activity correction is applied to the proton
   component when computing pH.  The vector must contain one entry per
   component (``NCOMP`` values).Components that should not contribute to the ionic 
   strength, should be assigned a charge of zero.

===================  =========================  =========================================
**Type:** int        **Range:** :math:`\ge 0`     **Length:** NCOMP
===================  =========================  =========================================

``CPA_MAXITER``
   Maximum number of Newton iterations for solving the adsorber surface
   potential :math:`\psi_{0,A}` (optional, defaults to 100). If the solver
   does not converge, a warning is emitted and the last iterate is used.

===================  =========================  =========================================
**Type:** int        **Range:** :math:`\ge 1`   **Length:** 1
===================  =========================  =========================================
