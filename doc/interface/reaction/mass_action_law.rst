.. _mass_action_law_config:

Mass Action Law
~~~~~~~~~~~~~~~

**Group /input/model/unit_XXX/reaction_phase/reaction_model_YYY - REACTION_MODEL = MASS_ACTION_LAW**

For information on model equations, refer to :ref:`mass_action_law_model`.

``MAL_KFWD``

   Forward rate constants for bulk volume reactions (available for external functions)
   
   ================  =========================  ===================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NREACT}`
   ================  =========================  ===================================
   
``MAL_KBWD``

   Backward rate constants for bulk volume reactions (available for external functions)
   
   ================  =========================  ===================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NREACT}`
   ================  =========================  ===================================

``MAL_STOICHIOMETRY``

   Stoichiometric matrix of bulk volume reactions as :math:`\texttt{NCOMP} \times \texttt{NREACT}` matrix in row-major storage
   
   ================  ========================================================
   **Type:** double  **Length:** :math:`\texttt{NCOMP} \cdot \texttt{NREACT}`
   ================  ========================================================

``MAL_EXPONENTS_FWD``

   Forward exponent matrix of bulk volume reactions as :math:`\texttt{NCOMP} \times \texttt{NREACT}` matrix in row-major storage (optional, calculated from :math:`\texttt{MAL_STOICHIOMETRY_BULK}` by default)
   
   ================  ========================================================
   **Type:** double  **Length:** :math:`\texttt{NCOMP} \cdot \texttt{NREACT}`
   ================  ========================================================

``MAL_EXPONENTS_BWD``

   Backward exponent matrix of bulk volume reactions as :math:`\texttt{NCOMP} \times \texttt{NREACT}` matrix in row-major storage (optional, calculated from :math:`\texttt{MAL_STOICHIOMETRY_BULK}` by default)
   
   ================  ========================================================
   **Type:** double  **Length:** :math:`\texttt{NCOMP} \cdot \texttt{NREACT}`
   ================  ========================================================
