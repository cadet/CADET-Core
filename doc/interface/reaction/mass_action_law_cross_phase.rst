.. _mass_action_law_cross_phase_config:

Mass Action Law Cross Phase
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Group /input/model/unit_XXX(/particle_type_YYY)/phase_reaction_ZZZ/ - REACTION_MODEL = MASS_ACTION_LAW_CROSS_PHASE**

For information on model equations, refer to :ref:`mass_action_law_model_cross_phase`.

   
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

Examples
--------

Cross-phase reaction that consumes one bulk and one liquid component and produces a solid state (one reaction)::

   input.model.unit_000.reaction_cross_phase_000.NREAC_CROSS_PHASE = 1
   input.model.unit_000.reaction_cross_phase_000.type = MASS_ACTION_LAW_CROSS_PHASE
   input.model.unit_000.reaction_cross_phase_000.reaction_model_000.MAL_KFWD_LIQUID = [1.0]
   input.model.unit_000.reaction_cross_phase_000.reaction_model_000.MAL_KFWD_SOLID = [1.0]
   input.model.unit_000.reaction_cross_phase_000.reaction_model_000.MAL_STOICHIOMETRY_LIQUID = [... length NCOMP*1 ...]
   input.model.unit_000.reaction_cross_phase_000.reaction_model_000.MAL_STOICHIOMETRY_SOLID = [... length NTOTALBND*1 ...]

Cross-Phase reaction in a particle::

   input.model.unit_000.particle_type_000.reaction_cross_phase_000.NREAC_CROSS_PHASE = 1
   input.model.unit_000.particle_type_000.reaction_cross_phase_000.type = MASS_ACTION_LAW_CROSS_PHASE
   input.model.unit_000.particle_type_000.reaction_cross_phase_000.reaction_model_000.MAL_KFWD_LIQUID = [1.0]
   input.model.unit_000.particle_type_000.reaction_cross_phase_000.reaction_model_000.MAL_KFWD_SOLID = [1.0]
   input.model.unit_000.particle_type_000.reaction_cross_phase_000.reaction_model_000.MAL_STOICHIOMETRY_LIQUID = [... length NCOMP*1 ...]
   input.model.unit_000.particle_type_000.reaction_cross_phase_000.reaction_model_000.MAL_STOICHIOMETRY_SOLID = [... length NTOTALBND*1 ...]



See also
--------

- Single-phase variant: :ref:`mass_action_law_config`
