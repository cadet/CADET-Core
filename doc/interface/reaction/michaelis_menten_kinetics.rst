.. _michaelis_menten_kinetics_config:

Michaelis Menten kinetics
~~~~~~~~~~~~~~~~~~~~~~~~~

**Group /input/model/unit_XXX/reaction - REACTION_MODEL = MICHAELIS_MENTEN**

For information on model equations, refer to :ref:`michaelis_menten_kinetics_model`.

``MM_STOICHIOMETRY_BULK``

   Stochiometric matrix :math:`S`.
   The substrate component :math:`c_S` is identified by the index of the first negative entry in the stoichiometry of the corresponding reaction.
   Input as reaction index major.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NREACT} \cdot \texttt{NCOMP}`
   ================  =============================  ========================================================

``MM_VMAX``

	Limiting rate :math:`\mu_{\mathrm{max},j}` at saturation.

   ================  =============================  ===================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NREACT}`
   ================  =============================  ===================================

``MM_KMM``

	Michaelis constant :math:`k_{\mathrm{MM},j}`.

   ================  =============================  ===================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NREACT}`
   ================  =============================  ===================================

``MM_KI_C``

	Inhibition constant for competitive inhibition :math:`K_{i,j,k}` w.r.t substrate :math:`i`, inhibitor :math:`j` and reaction :math:`j`. If :math:`k_{i,j,k} \leq 0`, the component does not inhibit the reaction.
	Input as reaction index major.

   ================  =============================  =============================================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NREACT} \cdot \texttt{NCOMP} \cdot \texttt{NCOMP}`
   ================  =============================  =============================================================================

``MM_KI_NC``

	Inhibition constant non competitive inhibition:math:`\tilde{K}_{i,j,k}` w.r.t substrate :math:`i`, inhibitor :math:`j` and reaction :math:`j`. If :math:`k_{i,j,k} \leq 0`, the component does not inhibit the reaction.
	Input as reaction index major.

   ================  =============================  =============================================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NREACT} \cdot \texttt{NCOMP} \cdot \texttt{NCOMP}`
   ================  =============================  =============================================================================

Example configuration
^^^^^^^^^^^^^^^^^^^^^
This example shows how to configure a Michaelis Menten reaction system in CADET-Python.
For a system with two components A and B, where A is the substrate and B is the product with 
* Michaelis constant of `km_a`,
* competitive inhibition constant of `ki_b` for B inhibiting A
* and a maximum rate of `vmax`.

.. code-block:: Python3

   #Configure the reaction system
model.root.input.model.unit_001.reaction_model = 'MICHAELIS_MENTEN'
        
# Km values 2D array [reaction][components]
model.root.input.model.unit_001.reaction_bulk.mm_kmm = [
    [km_a, 0.0] # A is substrate
]

# Competitive inhibition constants - 3D array [reaction][components][components]
model.root.input.model.unit_001.reaction_bulk.mm_ki_c = [
    [
        [0.0, ki_b], # Inhibition konstant for A (Product inhibtion B inhibits A)
        [0.0, 0.0],  # Inhibition konstant for B
    ]
]

# Uncompetitive inhibition constants - 3D array [reaction][components][components]
model.root.input.model.unit_001.reaction_bulk.mm_ki_uc = [
    [
        [0.0, 0.0], # Inhibition konstant for A 
        [0.0, 0.0], # Inhibition konstant for B
    ]
]

# Vmax values 1D array [reaction]
model.root.input.model.unit_001.reaction_bulk.mm_vmax = [vmax]

# Stoichiometry matrix 2D array [components][reaction]
model.root.input.model.unit_001.reaction_bulk.mm_stoichiometry_bulk = [
    [-1],
    [1] # A -> B
]
