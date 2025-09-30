.. _non_consistency_solver_parameters:

Nonlinear solver for consistent initialization
===============================================

Group /input/model/unit_XXX/discretization/consistency_solver - Nonlinear consistency solver paramters
------------------------------------------------------------------------------------------------------

The following specifications of the nonlinear solver for consistent initialization are available.
The default is a ``COMPOSITE`` solver that applies multiple solvers subsequently, where the default is ``ATRN_ERR`` followed by ``LEVMAR``.


``SOLVER_NAME``

Name of the solver.
Available solvers are

- ``LEVMAR``: Levenberg-Marquardt method (ported from the MATLAB implementation)

- ``ATRN_RES``: Adaptive trust-region Newton solver, implemetation of the NLEQ-RES algorithm described in :cite:`Deuflhard2011` (p. 131)

- ``ATRN_ERR``: Robust adaptive trust-region Newton solver, implemetation of the NLEQ-ERR algorithm described in :cite:`Deuflhard2011` (pp. 148)

- ``COMPOSITE``: Applies multiple solvers subsequently


  ================== =======================
   **Type:** string  **Length:** :math:`1`     
  ================== =======================

``SUBSOLVERS``

Vector with names of solvers for the ``COMPOSITE`` solver (only required for composite solver). See ``SOLVER_NAME`` for available solvers.

  ================== ==========================
   **Type:** string  **Length:** :math:`\gt 1`     
  ================== ==========================

``INIT_DAMPING``

Initial damping factor (default is :math:`0.01`), only required for ``LEVMAR``, ``ATRN_RES`` and ``ATRN_ERR``. See ``SOLVER_NAME`` for available solvers.

   ================  =============================  ==================================
   **Type:** double  **Range:** :math:`\ge 0`       **Length:** :math:`1`
   ================  =============================  ==================================

``MIN_DAMPING``

Minimal damping factor (default is :math:`0.0001`), only required for ``ATRN_RES`` and ``ATRN_ERR``. See ``SOLVER_NAME`` for available solvers.

   ================  =============================  ==================================
   **Type:** double  **Range:** :math:`\ge 0`       **Length:** :math:`1`
   ================  =============================  ==================================
