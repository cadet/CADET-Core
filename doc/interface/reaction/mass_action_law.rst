.. _mass_action_law_config:

Mass Action Law
~~~~~~~~~~~~~~~

**Group /input/model/unit_XXX(/particle_type_YYY)/phase_reaction_ZZZ/ - REACTION_MODEL = MASS_ACTION_LAW**

For information on model equations, refer to :ref:`mass_action_law_model`.

Notes
-----

- ``reaction_phase`` refers to one of the phase-specific groups listed in :ref:`FFReaction`, e.g., ``reaction_bulk``, ``reaction_solid``, or ``reaction_particle_YYY`` (for particle type ``YYY``).
- Each ``reaction_model_YYY`` is one instance of a reaction model and can contain multiple reactions (lengths denoted by ``NREACT`` below).
- Dimensions of matrices depend on the hosting phase of this model instance:

  - Bulk phase or particle liquid phase: ``NVAR = NCOMP``
   - Particle solid phase: ``NVAR = NTOTALBOUND`` (total number of bound states across all components)

``MAL_KFWD``

   Forward rate constants for reactions in this phase (available for external functions)
   
   ================  =========================  ===================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NREACT}`
   ================  =========================  ===================================
   
``MAL_KBWD``

   Backward rate constants for reactions in this phase (available for external functions)
   
   ================  =========================  ===================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NREACT}`
   ================  =========================  ===================================

``MAL_STOICHIOMETRY``

   Stoichiometric matrix as :math:`\texttt{NVAR} \times \texttt{NREACT}` matrix in row-major storage (phase-dependent ``NVAR``, see Notes)
   
   ================  ========================================================
   **Type:** double  **Length:** :math:`\texttt{NVAR} \cdot \texttt{NREACT}`
   ================  ========================================================

``MAL_EXPONENTS_FWD``

   Forward exponent matrix as :math:`\texttt{NVAR} \times \texttt{NREACT}` matrix in row-major storage (optional, calculated from :math:`\texttt{MAL_STOICHIOMETRY}` by default)
   
   ================  ========================================================
   **Type:** double  **Length:** :math:`\texttt{NVAR} \cdot \texttt{NREACT}`
   ================  ========================================================

``MAL_EXPONENTS_BWD``

   Backward exponent matrix as :math:`\texttt{NVAR} \times \texttt{NREACT}` matrix in row-major storage (optional, calculated from :math:`\texttt{MAL_STOICHIOMETRY}` by default)
   
   ================  ========================================================
   **Type:** double  **Length:** :math:`\texttt{NVAR} \cdot \texttt{NREACT}`
   ================  ========================================================

Examples
--------
.. code-block::

Cross-phase reaction that consumes one bulk and one liquid component and produces a solid state (one reaction)::

   input.model.unit_000.reaction_liquid_000.NREAC_LIQUID = 1
   input.model.unit_000.reaction_liquid_000.type = MASS_ACTION_LAW
   input.model.unit_000.reaction_liquid_000.reaction_model_000.MAL_KFWD = [1.0]
   input.model.unit_000.reaction_liquid_000.reaction_model_000.MAL_KBWD = [1.0]
   input.model.unit_000.reaction_liquid_000.reaction_model_000.MAL_STOICHIOMETRY = [... length NCOMP*1 ...]

Cross-Phase reaction in a pore::

   input.model.unit_000.particle_type_000.NREAC_PORE = 1
   input.model.unit_000.particle_type_000.reaction_pore_000.type = MASS_ACTION_LAW
   input.model.unit_000.particle_type_000.reaction_pore_000.reaction_model_000.MAL_KFWD = [1.0]
   input.model.unit_000.particle_type_000.reaction_pore_000.reaction_model_000.MAL_KBWD = [1.0]
   input.model.unit_000.particle_type_000.reaction_pore_000.reaction_model_000.MAL_STOICHIOMETRY = [... length NCOMP*1 ...]

See also
--------

- Cross-phase variant: :ref:`mass_action_law_cross_phase_config` (supports liquid/solid modifiers and bulk coupling)
