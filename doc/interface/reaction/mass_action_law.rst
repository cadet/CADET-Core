.. _mass_action_law_config:

Mass Action Law
~~~~~~~~~~~~~~~

**Group /input/model/unit_XXX/reaction - REACTION_MODEL = MASS_ACTION_LAW**

``MAL_KFWD_BULK``

   Forward rate constants for bulk volume reactions (available for external functions)
   
   ================  =========================  ===================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NREACT}`
   ================  =========================  ===================================
   
``MAL_KBWD_BULK``

   Backward rate constants for bulk volume reactions (available for external functions)
   
   ================  =========================  ===================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NREACT}`
   ================  =========================  ===================================
   
``MAL_KFWD_LIQUID``

   Forward rate constants for particle liquid phase reactions (available for external functions)
   
   ================  =========================  ===================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NREACT}`
   ================  =========================  ===================================
   
``MAL_KBWD_LIQUID``

   Backward rate constants for particle liquid phase reactions (available for external functions)
   
   ================  =========================  ===================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NREACT}`
   ================  =========================  ===================================
   
``MAL_KFWD_SOLID``

   Forward rate constants for particle solid phase reactions (available for external functions)
   
   ================  =========================  ===================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NREACT}`
   ================  =========================  ===================================
   
``MAL_KBWD_SOLID``

   Backward rate constants for particle solid phase reactions (available for external functions)
   
   ================  =========================  ===================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NREACT}`
   ================  =========================  ===================================
   
``MAL_STOICHIOMETRY_BULK``

   Stoichiometric matrix of bulk volume reactions as :math:`\texttt{NCOMP} \times \texttt{NREACT}` matrix in row-major storage
   
   ================  ========================================================
   **Type:** double  **Length:** :math:`\texttt{NCOMP} \cdot \texttt{NREACT}`
   ================  ========================================================
   
``MAL_EXPONENTS_BULK_FWD``

   Forward exponent matrix of bulk volume reactions as :math:`\texttt{NCOMP} \times \texttt{NREACT}` matrix in row-major storage (optional, calculated from :math:`\texttt{MAL_STOICHIOMETRY_BULK}` by default)
   
   ================  ========================================================
   **Type:** double  **Length:** :math:`\texttt{NCOMP} \cdot \texttt{NREACT}`
   ================  ========================================================
   
``MAL_EXPONENTS_BULK_BWD``

   Backward exponent matrix of bulk volume reactions as :math:`\texttt{NCOMP} \times \texttt{NREACT}` matrix in row-major storage (optional, calculated from :math:`\texttt{MAL_STOICHIOMETRY_BULK}` by default)
   
   ================  ========================================================
   **Type:** double  **Length:** :math:`\texttt{NCOMP} \cdot \texttt{NREACT}`
   ================  ========================================================
   
``MAL_STOICHIOMETRY_LIQUID``

   Stoichiometric matrix of particle liquid phase reactions as :math:`\texttt{NCOMP} \times \texttt{NREACT}` matrix in row-major storage
   
   ================  ========================================================
   **Type:** double  **Length:** :math:`\texttt{NCOMP} \cdot \texttt{NREACT}`
   ================  ========================================================
   
``MAL_EXPONENTS_LIQUID_FWD``

   Forward exponent matrix of particle liquid phase reactions as :math:`\texttt{NCOMP} \times \texttt{NREACT}` matrix in row-major storage (optional, calculated from :math:`\texttt{MAL_STOICHIOMETRY_LIQUID}` by default)
   
   ================  ========================================================
   **Type:** double  **Length:** :math:`\texttt{NCOMP} \cdot \texttt{NREACT}`
   ================  ========================================================
   
``MAL_EXPONENTS_LIQUID_BWD``

   Backward exponent matrix of particle liquid phase reactions as :math:`\texttt{NCOMP} \times \texttt{NREACT}` matrix in row-major storage (optional, calculated from :math:`\texttt{MAL_STOICHIOMETRY_LIQUID}` by default)
   
   ================  ========================================================
   **Type:** double  **Length:** :math:`\texttt{NCOMP} \cdot \texttt{NREACT}`
   ================  ========================================================
   
``MAL_EXPONENTS_LIQUID_FWD_MODSOLID``

   Forward solid phase modifier exponent matrix of particle liquid phase reactions as :math:`\texttt{NTOTALBND} \times \texttt{NREACT}` matrix in row-major storage (optional, defaults to all 0)
   
   ================  ============================================================
   **Type:** double  **Length:** :math:`\texttt{NTOTALBND} \cdot \texttt{NREACT}`
   ================  ============================================================
   
``MAL_EXPONENTS_LIQUID_BWD_MODSOLID``

   Backward solid phase modifier exponent matrix of particle liquid phase reactions as :math:`\texttt{NTOTALBND} \times \texttt{NREACT}` matrix in row-major storage (optional, defaults to all 0)
   
   ================  ============================================================
   **Type:** double  **Length:** :math:`\texttt{NTOTALBND} \cdot \texttt{NREACT}`
   ================  ============================================================
   
``MAL_STOICHIOMETRY_SOLID``

   Stoichiometric matrix of particle solid phase reactions as :math:`\texttt{NTOTALBND} \times \texttt{NREACT}` matrix in row-major storage
   
   ================  ============================================================
   **Type:** double  **Length:** :math:`\texttt{NTOTALBND} \cdot \texttt{NREACT}`
   ================  ============================================================
   
``MAL_EXPONENTS_SOLID_FWD``

   Forward exponent matrix of particle solid phase reactions as :math:`\texttt{NTOTALBND} \times \texttt{NREACT}` matrix in row-major storage (optional, calculated from :math:`\texttt{MAL_STOICHIOMETRY_SOLID}` by default)
   
   ================  ============================================================
   **Type:** double  **Length:** :math:`\texttt{NTOTALBND} \cdot \texttt{NREACT}`
   ================  ============================================================
   
``MAL_EXPONENTS_SOLID_BWD``

   Backward exponent matrix of particle solid phase reactions as :math:`\texttt{NTOTALBND} \times \texttt{NREACT}` matrix in row-major storage (optional, calculated from :math:`\texttt{MAL_STOICHIOMETRY_SOLID}` by default)
   
   ================  ============================================================
   **Type:** double  **Length:** :math:`\texttt{NTOTALBND} \cdot \texttt{NREACT}`
   ================  ============================================================
   
``MAL_EXPONENTS_SOLID_FWD_MODLIQUID``

   Forward liquid phase modifier exponent matrix of particle solid phase reactions as :math:`\texttt{NCOMP} \times \texttt{NREACT}` matrix in row-major storage (optional, defaults to all 0)
   
   ================  ========================================================
   **Type:** double  **Length:** :math:`\texttt{NCOMP} \cdot \texttt{NREACT}`
   ================  ========================================================
   
``MAL_EXPONENTS_SOLID_BWD_MODLIQUID``

   Backward liquid phase modifier exponent matrix of particle solid phase reactions as :math:`\texttt{NCOMP} \times \texttt{NREACT}` matrix in row-major storage (optional, defaults to all 0)
   
   ================  ========================================================
   **Type:** double  **Length:** :math:`\texttt{NCOMP} \cdot \texttt{NREACT}`
   ================  ========================================================
