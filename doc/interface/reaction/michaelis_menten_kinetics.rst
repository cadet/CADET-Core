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

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NREACT} \cdot \texttt{NCOMP} \cdot \texttt{NCOMP}`
   ================  =============================  ========================================================

``MM_KI_NC``

	Inhibition constant non competitive inhibition:math:`\tilde{K}_{i,j,k}` w.r.t substrate :math:`i`, inhibitor :math:`j` and reaction :math:`j`. If :math:`k_{i,j,k} \leq 0`, the component does not inhibit the reaction.
	Input as reaction index major.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NREACT} \cdot \texttt{NCOMP} \cdot \texttt{NCOMP}`
   ================  =============================  ========================================================
