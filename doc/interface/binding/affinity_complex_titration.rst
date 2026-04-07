.. _affinity_complex_titration_config:

Affinity Complex Titration
==========================

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = AFFINITY_COMPLEX_TITRATION**

For information on model equations, refer to :ref:`affinity_complex_titration`. 

General remarks
---------------

* The **first component** is the ACT modulator.
* This first component must be **non-binding**.
* The ACT implementation currently supports **at most one bound state per component**.
* The first component can be specified either as a negative logarithmic quantity
  (:math:`\mathrm{pH}`, :math:`\mathrm{pNa}`, ...) or as a raw ion concentration.

``IS_KINETIC``
   Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 = quasi-stationary.
   If a single value is given, the mode is set for all bound states. Otherwise, the adsorption
   mode is set for each bound state separately.

===================  =========================  =========================================
**Type:** int        **Range:** {0,1}           **Length:** 1/NTOTALBND
===================  =========================  =========================================

``ACT_USE_ION_CONC``
   Selects how the first liquid-phase component is interpreted.

   * ``False``: the first component is already on the logarithmic axis, for example
     :math:`\mathrm{pH}` or :math:`\mathrm{pNa}`.
     Use ``ACT_PKAA`` and ``ACT_PKAG``.
   * ``True``: the first component is the raw ion concentration :math:`c_{\mathrm{ion}}`.
     Use ``ACT_CMID_A`` and ``ACT_CMID_G``.

   Default is ``False``.

   The two conventions are equivalent when

   .. math::

      \mathrm{p}K_{a,A,i} = -\log_{10}(c_{\mathrm{mid},A,i}), \qquad
      \mathrm{p}K_{a,G,i} = -\log_{10}(c_{\mathrm{mid},G,i}).

===================  =========================  =========================================
**Type:** bool       **Range:** {False, True}   **Length:** 1
===================  =========================  =========================================

Core kinetic parameters
-----------------------

``ACT_KA``
   Adsorption rate constants.

**Unit:** :math:`m_{MP}^3\,mol^{-1}\,s^{-1}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================

``ACT_KD``
   Desorption rate constants.

**Unit:** :math:`s^{-1}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================

``ACT_QMAX``
   Maximum binding capacities before modulation by the ACT capacity gate.

**Unit:** :math:`mol\,m_{SP}^{-3}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`> 0`     **Length:** NCOMP
===================  =========================  =========================================

``ACT_ETAA``
   Hill-type coefficients controlling how strongly the apparent binding capacity changes with the modulator.

**Unit:** :math:`1`

===================  =============================  =========================================
**Type:** double     **Range:** :math:`\mathbb{R}`   **Length:** NCOMP
===================  =============================  =========================================

``ACT_ETAG``
   Hill-type coefficients controlling how strongly the effective adsorption strength changes with the modulator.

**Unit:** :math:`1`

===================  =============================  =========================================
**Type:** double     **Range:** :math:`\mathbb{R}`   **Length:** NCOMP
===================  =============================  =========================================

Midpoint parameters for ``ACT_USE_ION_CONC = False``
----------------------------------------------------

Use these parameters only when the first component is given as :math:`\mathrm{pIon}`
(for example :math:`\mathrm{pH}`).

``ACT_PKAA``
   Midpoint of the capacity transition on the logarithmic ion axis.

**Unit:** :math:`1`

===================  ==============================  =========================================
**Type:** double     **Range:** :math:`\mathbb{R}`   **Length:** NCOMP
===================  ==============================  =========================================

``ACT_PKAG``
   Midpoint of the equilibrium / adsorption-strength transition on the logarithmic ion axis.

**Unit:** :math:`1`

===================  ==============================  =========================================
**Type:** double     **Range:** :math:`\mathbb{R}`   **Length:** NCOMP
===================  ==============================  =========================================

Midpoint parameters for ``ACT_USE_ION_CONC = True``
---------------------------------------------------

Use these parameters only when the first component is given as a raw ion concentration.
CADET internally converts them to the same negative logarithmic axis used by the ``ACT_PKAA`` / ``ACT_PKAG`` form.

``ACT_CMID_A``
   Midpoint ion concentration for the capacity transition.
   Must be non-negative and should use the same concentration unit as the first liquid-phase component.

**Recommended unit:** :math:`mol\,m_{MP}^{-3}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================

``ACT_CMID_G``
   Midpoint ion concentration for the equilibrium / adsorption-strength transition.
   Must be non-negative and should use the same concentration unit as the first liquid-phase component.

**Recommended unit:** :math:`mol\,m_{MP}^{-3}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================

Practical switching rule
------------------------

To convert an existing ACT setup from ``ACT_USE_ION_CONC = False`` to ``ACT_USE_ION_CONC = True`` while keeping
exactly the same model response:

1. Replace the first component value :math:`\mathrm{pIon}` by :math:`c_{\mathrm{ion}} = 10^{-\mathrm{pIon}}`.
2. Replace ``ACT_PKAA`` by ``ACT_CMID_A = 10^{-\mathrm{ACT\_PKAA}}``.
3. Replace ``ACT_PKAG`` by ``ACT_CMID_G = 10^{-\mathrm{ACT\_PKAG}}``.
4. Keep ``ACT_KA``, ``ACT_KD``, ``ACT_QMAX``, ``ACT_ETAA``, and ``ACT_ETAG`` unchanged.
