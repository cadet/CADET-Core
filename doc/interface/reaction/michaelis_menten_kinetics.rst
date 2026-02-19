.. _michaelis_menten_kinetics_config:

Michaelis Menten kinetics
~~~~~~~~~~~~~~~~~~~~~~~~~

**Group /input/model/unit_XXX(/particle_type_YYY)/phase_reaction_ZZZ/ - TYPE = MICHAELIS_MENTEN**

For information on model equations, refer to :ref:`michaelis_menten_kinetics_model`.

Notes
-----

- ``phase_reaction_ZZZ`` refers to one of the phase-specific raction groups listed in :ref:`FFReaction`, e.g., ``liquid_reaction_ZZZ``, ``solid_reaction_ZZZ``.
- Some dimensions below depend on the hosting phase of this model instance:
    - Bulk phase or particle liquid phase: ``NVAR = NCOMP``
    - Particle solid phase: ``NVAR = NTOTALBOUND`` (total number of bound states across all components)

``MM_STOICHIOMETRY_BULK``
    Stoichiometric matrix :math:`S`.
    This matrix defines the quantitative relationships between reactants and products for each reaction in the system.
    Each entry :math:`S_{i,j}` specifies the stoichiometric coefficient for component :math:`i` in reaction :math:`j`.
    Negative values indicate consumption (substrate), while positive values indicate production (products).
    Input as reaction index major.

    **Unit:** None

    ================  =============================  ========================================================
    **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NREACT} \cdot \texttt{NVAR}`
    ================  =============================  ========================================================

``MM_VMAX``
    Maximum reaction rate :math:`v_{\mathrm{max},j}` at substrate saturation for reaction :math:`j`.
    This parameter defines the upper limit of the reaction rate when the substrate concentration is sufficiently high such that the enzyme is saturated.

    **Unit:** :math:`~mol^{-1}~m^{-3}~s^{-1}`

    ================  =============================  ===================================
    **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NREACT}`
    ================  =============================  ===================================

``MM_KM``
    Michaelis constant :math:`K_{\mathrm{M}_{i,j}}` for reaction :math:`j` and substrate :math:`i`.
    This constant represents the substrate concentration at which the reaction rate is half of its maximum value.

    **Unit:** :math:`~mol^{-1}~m^{-3}`

    ================  =============================  ===================================
    **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NREACT}`
    ================  =============================  ===================================

``MM_KI_C``
    Inhibition constant for competitive inhibition :math:`K^{c}_{I_{k}}`.
    The index :math:`k` corresponds to the inhibitors acting on substrate :math:`c_{i,j}` in reaction :math:`j`, i.e. :math:`k = (j,i,k)`, where :math:`k` is the index of the inhibitor.
    If :math:`K^{c}_{I_{k}} > 0`, the component inhibits the reaction.
    Input as reaction index major.

    **Unit:** :math:`mol^{-1}~m^{-3}`

    ================  =============================  =============================================================================
    **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NREACT} \cdot \texttt{NVAR} \cdot \texttt{NVAR}`
    ================  =============================  =============================================================================

``MM_KI_UC``
    Inhibition constant for uncompetitive inhibition :math:`K^{uc}_{I_{k}}`.
    The index :math:`k` corresponds to the inhibitors acting on substrate :math:`c_{i,j}` in reaction :math:`j`, i.e. :math:`k = (j,i,k)`, where :math:`k` is the index of the inhibitor.
    Input as reaction index major.

    **Unit:** :math:`mol^{-1}~m^{-3}`

    ================  =============================  =============================================================================
    **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NREACT} \cdot \texttt{NVAR} \cdot \texttt{NVAR}`
    ================  =============================  =============================================================================

CADET Python Interface Example
------------------------------
This example shows the configuration of one Michaelis-Menten reaction system in bulk liquid phase.
The system has two components A and B, where A is the substrate and B is the product.
In addition to that the model includes:

* A Michaelis constant ``KM_a``,
* competitive inhibition constant of ``KI_b_a`` for B inhibiting A,
* and a maximum rate of ``vmax``.

.. code-block:: python3

    # Configure the reaction system
    input.model.unit_000.NREAC_LIQUID = 1
    input.model.unit_000.reaction_liquid_000.type = 'MICHAELIS_MENTEN'

    # Km values 2D array [reaction][components]
    input.model.unit_000.reaction_liquid_000.mm_km = [
        [KM_a, 0.0]  # A is substrate
    ]

    # Competitive inhibition constants - 3D array [reaction][components][components]
    input.model.unit_000.reaction_liquid_000.mm_ki_c = [
        [
            [0.0, KI_b_a],  # Inhibition constant for A (Product inhibition B inhibits A)
            [0.0, 0.0],     # Inhibition constant for B (not active)
        ]
    ]

    # Uncompetitive inhibition constants - 3D array [reaction][components][components]
    input.model.unit_000.reaction_liquid_000.mm_ki_uc = [
        [
            [0.0, 0.0],  # Inhibition constant for A (not active)
            [0.0, 0.0],  # Inhibition constant for B (not active)
        ]
    ]

    # Vmax values 1D array [reaction]
    input.model.unit_000.reaction_liquid_000.mm_vmax = [vmax]

    # Stoichiometry matrix 2D array [components][reaction]
    input.model.unit_000.reaction_liquid_000.mm_stoichiometry = [
        [-1],
        [1]  # A -> B
    ]
