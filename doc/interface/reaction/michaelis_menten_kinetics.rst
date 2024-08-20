.. _michaelis_menten_kinetics_config:

Michaelis Menten kinetics
~~~~~~~~~~~~~~~~~~~~~~~~~

**Group /input/model/unit_XXX/reaction - REACTION_MODEL = MICHAELIS_MENTEN**

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

``MM_KI``

	Inhibition constant :math:`k_{\mathrm{I},j,i}` w.r.t component :math:`i` and reaction :math:`j`. If :math:`k_{\mathrm{I},j,i} <= 0`, the component does not inhibit the reaction.
	Input as reaction index major.
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NREACT} \cdot \texttt{NCOMP}`
   ================  =============================  ========================================================