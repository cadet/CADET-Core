.. _michaelis_menten_kinetics_config:

Michaelis Menten kinetics
~~~~~~~~~~~~~~~~~~~~~~~~~

**Group /input/model/unit_XXX/reaction - REACTION_MODEL = MICHAELIS_MENTEN**

For information on model equations, refer to :ref:`_michaelis_menten_kinetics_model`.

``MM_STOICHIOMETRY_BULK``

    Stoichiometric matrix :math:`S`.
    This matrix defines the quantitative relationships between reactants and products for each reaction in the system.
    Each entry :math:`S_{i,j}` specifies the stoichiometric coefficient for component :math:`i` in reaction :math:`j`.
    Negative values indicate consumption (reactants), while positive values indicate production (products).
   Input as reaction index major.
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NREACT} \cdot \texttt{NCOMP}`
   ================  =============================  ========================================================
   
``MM_VMAX``

    Maximum reaction rate :math:`v_{\mathrm{max},j}` at substrate saturation for reaction :math:`j`.
    This parameter defines the upper limit of the reaction rate when the substrate concentration is sufficiently high such that the enzyme is saturated.
   
    **Unit:** :math:`~mol^{-1}~m^{-3}~s^{-1}`

   ================  =============================  ===================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NREACT}`
   ================  =============================  ===================================

``MM_KM``

    Michaelis constant :math:`K_{\mathrm{M},j}` for reaction :math:`j`.
    This constant represents the substrate concentration at which the reaction rate is half of its maximum value (:math:`v_{\mathrm{max},j}`).
   
    **Unit:** :math:`~mol^{-1}~m^{-3}`

   ================  =============================  ===================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NREACT}`
   ================  =============================  ===================================

``MM_KI_C``

	Inhibition constant for competitive inhibition :math:`K_{I_{i,j,k}}` w.r.t substrate :math:`i`, inhibitor :math:`j` and reaction :math:`j`. If :math:`K_{I_{i,j,k}} \leq 0`, the component does not inhibit the reaction.
	Input as reaction index major.
   
    **Unit:** :math:`mol^{-1}~m^{-3}`

   ================  =============================  =============================================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NREACT} \cdot \texttt{NCOMP} \cdot \texttt{NCOMP}`
   ================  =============================  =============================================================================

``MM_KI_NC``

	Inhibition constant non competitive inhibition :math:`\tilde{K}_{I_{i,j,k}}` w.r.t substrate :math:`i`, inhibitor :math:`j` and reaction :math:`j`. If :math:`\tilde{K}_{I_{i,j,k}} \leq 0`, the component does not inhibit the reaction.
	Input as reaction index major.

    **Unit:** :math:`mol^{-1}~m^{-3}`

   ================  =============================  =============================================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NREACT} \cdot \texttt{NCOMP} \cdot \texttt{NCOMP}`
   ================  =============================  =============================================================================

Example configuration
^^^^^^^^^^^^^^^^^^^^^
This example shows the configuration of a Michaelis-Menten reaction system in CADET-Python.
The system has two components A and B, where A is the substrate and B is the product. 
In addition to that the model includes:

* a Michaelis constant ``km_a``,
* competitive inhibition constant of ``ki_b`` for B inhibiting A,
* and a maximum rate of ``vmax``

.. code-block:: Python3

   #Configure the reaction system
    model.root.input.model.unit_001.reaction_model = 'MICHAELIS_MENTEN'
            
    # Km values 2D array [reaction][components]
    model.root.input.model.unit_001.reaction_bulk.mm_km = [
        [km_a, 0.0] # A is substrate
    ]

    # Competitive inhibition constants - 3D array [reaction][components][components]
    model.root.input.model.unit_001.reaction_bulk.mm_ki_c = [
        [
            [0.0, ki_b], # Inhibition konstant for A (Product inhibtion B inhibits A)
            [0.0, 0.0],  # Inhibition konstant for B (not active)
        ]
    ]

    # Uncompetitive inhibition constants - 3D array [reaction][components][components]
    model.root.input.model.unit_001.reaction_bulk.mm_ki_uc = [
        [
            [0.0, 0.0], # Inhibition konstant for A (not active)
            [0.0, 0.0], # Inhibition konstant for B (not active)
        ]
    ]

    # Vmax values 1D array [reaction]
    model.root.input.model.unit_001.reaction_bulk.mm_vmax = [vmax]

    # Stoichiometry matrix 2D array [components][reaction]
    model.root.input.model.unit_001.reaction_bulk.mm_stoichiometry_bulk = [
        [-1],
        [1] # A -> B
    ]
