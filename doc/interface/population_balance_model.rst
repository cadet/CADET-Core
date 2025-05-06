.. _pbm_config:

Crystallization / Precipitation models
======================================

Crystallization / Precipitation in CADET can be modeled by :ref:`pbm_model`, :ref:`aggregation_model`, :ref:`fragmentation_model`, as well as their combinations.
The configuration is described for all three modules in the following.
For more information on the model equations, please refer to :ref:`FFCrystallization`.

Note that any of these models can be used in any unit operation that supports reactions.

.. attention::

   By setting the ``NCOMP`` field, you specify the number of FV cells for the internal coordinate:
   If the PBM is considered, we have two components that account for the solute :math:`c` and solubility :math:`c_\text{eq}`.
   Note that the first component must be solute :math:`c` and the last component must be the solubility :math:`c_\text{eq}`.
   Additionally, we discretize the particle size domain (internal coordinate) by some customary FV methods, giving us a finite set of particle sizes under consideration :math:`\{x_1, \dots, x_{N_x}\}`.
   Every particle size considered is treated as an individual component of the unit operation and the field ``NCOMP`` of that unit operation in which the crystallization happens, must be specified accordingly as :math:`N_x` or :math:`N_x + 2` if the model includes the PBM.

Example code for configuring the crystallization models is available in `CADET-Verification <https://github.com/cadet/CADET-Verification/>`_ .

Group /input/model/unit_XXX
---------------------------

``NCOMP``

   Number of components, which is defined by the number of considered particle sizes (i.e. discretization with FV cells) plus, in the case that the PBM is employed, 2 (for solute and solubility).

   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 3`  **Length:** 1
   =============  =========================  =============

``REACTION_MODEL``

   The crystallization code is implemented as a reaction module, which is why crystallization needs to be specified here

   ================  ==============================================  =============
   **Type:** String  **Range:** :math:`\texttt{ CRYSTALLIZATION \}`  **Length:** 1
   ================  ==============================================  =============

Group /input/model/unit_XXX/reaction_bulk - REACTION_MODEL = CRYSTALLIZATION - UNIT_TYPE = CSTR
-----------------------------------------------------------------------------------------------

*The following parameters need to be specified under Group /input/model/unit_XXX/reaction/reaction_bulk/ for CSTR units, and Group /input/model/unit_XXX/reaction/reaction/ for transport units like the LRM.*

``CRY_MODE``

   Crystallization mode, which determines the exact model equation to be employed.

  1. Pure PBM, as described in :ref:`pbm_model`.
  2. Pure aggregation, as described in
  3. Combined PBM and aggregation.
  4. Pure fragmentation, as described in
  5. Combined PBM and fragmentation.
  6. Combined aggregation and fragmentation.
  7. Combined PBM, aggregation, and fragmentation.

   ================  =========================  =====================================
   **Type:** double  **Range:** :math:`\geq 1`   **Length:** :math:`\mathrm{N_x} + 1`
   ================  =========================  =====================================

``CRY_BINS``

   Coordinates of the cell faces, e.g. equidistant or logarithmic discretization of the internal coordinate :math:`x \in [x_c, x_\text{end}]`, including the end points.

   ================  =========================  =====================================
   **Type:** double  **Range:** :math:`\geq 1`   **Length:** :math:`\mathrm{N_x} + 1`
   ================  =========================  =====================================

Population Mass Balance input parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``CRY_NUCLEI_MASS_DENSITY``

   Nulcei mass density

   ================  =========================  =============
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** 1
   ================  =========================  =============

``CRY_VOL_SHAPE_FACTOR``

   Volumetric shape factor of the particles

   ================  =========================  =============
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** 1
   ================  =========================  =============

``CRY_PRIMARY_NUCLEATION_RATE``

   Primary nucleation rate constant :math:`k_p`

   ================  =========================  =============
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** 1
   ================  =========================  =============

``CRY_SECONDARY_NUCLEATION_RATE``

   Secondary nucleation rate :math:`k_b`

   ================  =========================  =============
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** 1
   ================  =========================  =============

``CRY_GROWTH_RATE_CONSTANT``

   Growth rate constant :math:`k_g`

   ================  =========================  =============
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** 1
   ================  =========================  =============

``CRY_GROWTH_CONSTANT``

   Growth constant :math:`\gamma`

   ================  =========================  =============
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** 1
   ================  =========================  =============

``CRY_A``

   Defines constant :math:`a` used to determine the growth rate

   ================  =========================  =============
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** 1
   ================  =========================  =============

``CRY_G``

   Defines constant :math:`g` used to determine the growth rate

   ================  =========================  =============
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** 1
   ================  =========================  =============

``CRY_P``

   Defines constant :math:`p`  used to determine the growth rate

   ================  =========================  =============
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** 1
   ================  =========================  =============

``CRY_GROWTH_DISPERSION_RATE``

   Growth dispersion rate :math:`D_g`

   ================  =========================  =============
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** 1
   ================  =========================  =============

``CRY_U``

   Defines constant :math:`u` used to determine the primary nucleation

   ================  =========================  =============
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** 1
   ================  =========================  =============

``CRY_B``

   Defines constant :math:`b` used to determine the secondary nucleation

   ================  =========================  =============
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** 1
   ================  =========================  =============

``CRY_K``

   Defines constant :math:`k` used to determine the secondary nucleation

   ================  =========================  =============
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** 1
   ================  =========================  =============

``CRY_GROWTH_SCHEME_ORDER``

   Defines the growth flux FV reconstruction scheme. It can only be

   - :math:`1`: upwind scheme
   - :math:`2`: HR Koren scheme
   - :math:`3`: WENO23 scheme
   - :math:`4`: WENO35 scheme.

   We recommend using the HR Koren scheme, which showed to be the most performant in our benchmarks.

   =============  ================================  =============
   **Type:** int  **Range:** :math:`[1, \dots, 4]`  **Length:** 1
   =============  ================================  =============

Aggregation input parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``CRY_AGGREGATION_INDEX``

   Defines the aggregation kernel. It can only be

   - :math:`0`: constant kernel
   - :math:`1`: Brownian kernel
   - :math:`2`: Smoluchowski kernel
   - :math:`3`: Golovin kernel
   - :math:`4`: differential force kernel

   =============  ================================  =============
   **Type:** int  **Range:** :math:`[0, \dots, 4]`  **Length:** 1
   =============  ================================  =============

``CRY_AGGREGATION_RATE_CONSTANT``

   Aggregation rate constant :math:`\beta_0`

   ================  ======================  =============
   **Type:** double  **Range:** :math:`> 0`  **Length:** 1
   ================  ======================  =============


Fragmentation input parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``CRY_FRAGMENTATION_RATE_CONSTANT``

   Fragmentation rate constant :math:`S_0`

   ================  ======================  =============
   **Type:** double  **Range:** :math:`> 0`  **Length:** 1
   ================  ======================  =============

``CRY_FRAGMENTATION_KERNEL_GAMMA``

   Fragmentation kernel coefficient :math:`\gamma`

   ================  ========================  =============
   **Type:** double  **Range:** :math:`> 1.0`  **Length:** 1
   ================  ========================  =============

``CRY_FRAGMENTATION_SELECTION_FUNCTION_ALPHA``

   Fragmentation selection function coefficient :math:`\alpha`

   ================  ======================  =============
   **Type:** double  **Range:** :math:`> 0`  **Length:** 1
   ================  ======================  =============
