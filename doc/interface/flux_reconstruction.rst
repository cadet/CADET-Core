.. _flux_reconstruction_methods:

Flux reconstruction methods
===========================

Group /input/model/unit_XXX/discretization/reconstruction = WENO
-----------------------------------------------------------------

``BOUNDARY_MODEL``

    Boundary model type:
    0. Lower WENO order (stable)
    1. Zero weights (unstable for small :math:`D_{\mathrm{ax}}`)
    2. Zero weights for :math:`p \neq 0` (less stable)
    
    =============  ==============================  =============
    **Type:** int  **Range:** :math:`\{0, 1, 2\}`  **Length:** 1
    =============  ==============================  =============

``WENO_EPS``

    WENO :math:`\varepsilon`
    
    ================  ==================================  =============
    **Type:** double  **Range:** :math:`\mathbb{R}^{>0}`  **Length:** 1
    ================  ==================================  =============

``WENO_ORDER``

   WENO order, also called WENO :math:`k`:

   1. Standard upwind scheme (order 1)
   2. WENO 2 (order 3)
   3. WENO 3 (order 5)
   
   =============  ==============================  =============
   **Type:** int  **Range:** :math:`\{1, 2, 3\}`  **Length:** 1
   =============  ==============================  =============


Group /input/model/unit_XXX/discretization/reconstruction = KOREN
-----------------------------------------------------------------

The Koren scheme implemented in CADET intrinsically uses a van Leer flux limiter. It can reach a maximum order of 2 depending on the smoothness of the solution. The
BOUNDARY_MODEL is intrinsically set to 0 (stable).

``KOREN_EPS``

   Set :math:`\varepsilon` in the van Leer flux limiter

   ================  =========================  =============
   **Type:** double  **Range:** :math:`\qe 0\`  **Length:** 1
   ================  =========================  =============