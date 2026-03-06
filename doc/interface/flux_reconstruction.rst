.. _flux_reconstruction_methods:

Flux reconstruction methods
===========================


Group /input/model/unit_XXX/discretization
------------------------------------------

``GRID_FACES``

    Optional, required only for non-equidistant grids.
    An array specifying the coordinates of all grid faces, which are used for reconstruction.
    The array length must equal the number of grid faces, i.e., number of grid cells plus one, and the column boundaries must be included.

    .. note::
        - large ratios of adjacent cells sizes **may** reduce reconstruction accuracy. Recommendation :math:`r = max(Δx_{i+1}/Δx_i, Δx_i/Δx_{i+1}) < 3`.
        - Avoid large ratios :math:`Δx_{max} / Δx_{min}`, which **might** cause stiff ODE systems and thus slow, unstable time integration.

   ================  =========================  ============================
   **Type:** double  **Range:** :math:`[0, L]`  **Length:** :math:`NCOL + 1`
   ================  =========================  ============================


Group /input/model/unit_XXX/discretization/weno - RECONSTRUCTION = WENO
-----------------------------------------------------------------------

``BOUNDARY_MODEL``

    Boundary model type:
    0. Lower WENO order (stable)
    1. Zero weights (unstable for small :math:`D_{\mathrm{ax}}`)
    2. Zero weights for :math:`p \neq 0` (less stable)

    =============  ==============================  =============
    **Type:** int  **Range:** :math:`\{0, 1, 2\}`  **Length:** 1
    =============  ==============================  =============

``WENO_EPS``

    WENO :math:`\varepsilon`, a small regularization parameter added to the nonlinear weight denominator to prevent division by zero and improve numerical condition for small smoothness indicators.
    A default value :math:`1e-10` often works well.

    ================  =====================  =============
    **Type:** double  **Range:** :math:`>0`  **Length:** 1
    ================  =====================  =============

``WENO_ORDER``

   WENO order, also called WENO :math:`k`:

   1. Standard upwind scheme (order 1)
   2. WENO 2 (order 3, 2nd order at boundaries)
   3. WENO 3 (order 5, 2nd order at boundaries)

   =============  ==============================  =============
   **Type:** int  **Range:** :math:`\{1, 2, 3\}`  **Length:** 1
   =============  ==============================  =============


Group /input/model/unit_XXX/discretization/koren - RECONSTRUCTION = KOREN
-------------------------------------------------------------------------

The Koren scheme implemented in CADET intrinsically uses a van Leer flux limiter and has a theoretical convergence order between 1 and 2 depending on the smoothness of the solution.

``KOREN_EPS``

   Sets :math:`\varepsilon` in the van Leer flux limiter, a small numerical parameter added to the slope-ratio denominator to prevent division by zero and improve numerical condition when gradients are very small.
   A default value :math:`1e-10` often works well.

   ================  ======================  =============
   **Type:** double  **Range:** :math:`> 0`  **Length:** 1
   ================  ======================  =============
