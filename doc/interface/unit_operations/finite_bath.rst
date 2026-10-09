.. _finite_bath_config:

Finite bath
===========

Group /input/model/unit_XXX - UNIT_TYPE = FINITE_BATH
------------------------------------------------------

For information on model equations, refer to :ref:`finite_bath_model`.

The finite bath and the :ref:`cstr_config` share the same homogeneous bulk model and are selected by their particle types.
A unit operation that specifies :math:`\texttt{FINITE_BATH}` or :math:`\texttt{CSTR}` is configured as

- a CSTR, if all of its particle types are in rapid equilibrium with the bulk liquid, that is, :math:`\texttt{HAS_FILM_DIFFUSION} = 0`,
- a finite bath, if all of its particle types have a film diffusion resistance, that is, :math:`\texttt{HAS_FILM_DIFFUSION} = 1`.

Mixing both kinds of particle types in one unit operation is an error.
If :math:`\texttt{HAS_FILM_DIFFUSION}` is not given, it defaults to the unit type that was specified, that is, to 0 for :math:`\texttt{CSTR}` and to 1 for :math:`\texttt{FINITE_BATH}`.
Without particle types, the unit operation is configured as a CSTR if :math:`\texttt{INIT_LIQUID_VOLUME}` or :math:`\texttt{CONST_SOLID_VOLUME}` is given, and as a finite bath otherwise.

``UNIT_TYPE``

   Specifies the type of unit operation model

   ================  ========================================  =============
   **Type:** string  **Range:** :math:`\texttt{FINITE_BATH}`    **Length:** 1
   ================  ========================================  =============

``NCOMP``

   Number of chemical components

   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``NPARTYPE``

   Number of particle types (optional, defaults to 0). Each particle type is configured in its own group ``particle_type_XXX``, see :ref:`particle_model_config`.

   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 0`  **Length:** 1
   =============  =========================  =============

``LIQUID_VOLUME``

   Volume of the bulk liquid, which is constant. The flow rates into and out of the unit operation therefore have to be equal.

   **Unit:** :math:`\mathrm{m}^{3}`

   ================  =====================  =============
   **Type:** double  **Range:** :math:`>0`  **Length:** 1
   ================  =====================  =============

``BULK_POROSITY``

   Ratio of the bulk liquid volume to the total volume, :math:`\varepsilon_b = V^{\ell} / (V^{\ell} + V^{p})`, where :math:`V^{p}` is the total volume of the particles. Required if :math:`\texttt{NPARTYPE} \geq 1`, otherwise optional and defaulting to 1.

   ================  ========================  =============
   **Type:** double  **Range:** :math:`(0,1]`  **Length:** 1
   ================  ========================  =============

``PAR_TYPE_VOLFRAC``

   Volume fractions of the particle types, have to sum to 1. Required if :math:`\texttt{NPARTYPE} > 1`.

   ================  ========================  =====================================
   **Type:** double  **Range:** :math:`[0,1]`  **Length:** :math:`\texttt{NPARTYPE}`
   ================  ========================  =====================================

``INIT_C``

   Initial concentrations for each component in the bulk liquid phase

   **Unit:** :math:`\mathrm{mol}\,\mathrm{m}_{\mathrm{IV}}^{-3}`

   ================  =========================  ==================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NCOMP}`
   ================  =========================  ==================================

``INIT_STATE``

   Full state vector for initialization (optional, :math:`\texttt{INIT_C}`, :math:`\texttt{INIT_CP}` and :math:`\texttt{INIT_CS}` will be ignored; if length is :math:`2\texttt{NDOF}`, then the second half is used for time derivatives).
   The ordering of the state vector is defined in :ref:`UnitOperationStateOrdering`.

   **Unit:** :math:`various`

   ================  =============================  ====================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NDOF} / 2\texttt{NDOF}`
   ================  =============================  ====================================================

``NREAC_LIQUID``

   Number of bulk liquid phase reaction models (optional, only if liquid reactions are present).

   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 0`  **Length:** 1
   =============  =========================  =============


Group /input/model/unit_XXX/particle_type_XXX
----------------------------------------------

One group per particle type, holding the particle geometry, the film and pore transport, the binding model and the particle reactions, see :ref:`particle_model_config`.
A particle type of a finite bath has to specify :math:`\texttt{HAS_FILM_DIFFUSION} = 1`.
Since the bulk liquid of a well mixed vessel is not flowing, a velocity dependent :math:`\texttt{FILM_DIFFUSION_DEP}` is not supported.


Group /input/model/unit_XXX/discretization
-------------------------------------------

``USE_ANALYTIC_JACOBIAN``

   Determines whether the analytically computed Jacobian matrix (faster) is used (value is 1) instead of a Jacobian generated by algorithmic differentiation (slower, value is 0)

   =============  ===========================  =============
   **Type:** int  **Range:** :math:`\{0, 1\}`  **Length:** 1
   =============  ===========================  =============

``LINEAR_SOLVER``

   Sparse linear solver that is used for the linear systems of the unit operation (optional, defaults to :math:`\texttt{SparseLU}`)

   ================  ===================================================  =============
   **Type:** string  **Range:** :math:`\texttt{SparseLU}`, ...             **Length:** 1
   ================  ===================================================  =============

The group additionally holds the settings of the nonlinear solver used for consistent initialization, see :ref:`non_consistency_solver_parameters`.
The spatial discretization of a general rate particle type is given in its own group, see :ref:`particle_model_config`.


Group /input/model/unit_XXX/liquid_reaction_YYY
------------------------------------------------
Each bulk liquid phase reaction is specified in another subgroup `liquid_reaction_YYY`, see :ref:`FFReaction`.
Solid phase and cross-phase reactions are placed in the group of the respective particle type.
