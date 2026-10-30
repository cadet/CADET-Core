.. _cstr_config:

Continuous stirred tank reactor model
=====================================

Group /input/model/unit_XXX - UNIT_TYPE = CSTR
----------------------------------------------


For information on model equations, refer to :ref:`cstr_model`.

``UNIT_TYPE``

   Specifies the type of unit operation model

   ================  ================================  =============
   **Type:** string  **Range:** :math:`\texttt{CSTR}`  **Length:** 1
   ================  ================================  =============

``NCOMP``

   Number of chemical components in the chromatographic medium

   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``NPARTYPE``

   Number of particle types (optional, defaults to 0). Each particle type is configured in its own group ``particle_type_XXX``, see below.

   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 0`  **Length:** 1
   =============  =========================  =============

``INIT_C``

   Initial concentrations for each component in the bulk liquid phase

   **Unit:** :math:`\mathrm{mol}\,\mathrm{m}_{\mathrm{IV}}^{-3}`

   ================  =========================  ==================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NCOMP}`
   ================  =========================  ==================================

The volumes are given by one of two mutually exclusive parameterizations.
Specifying fields of both is an error.

1. Variable volume: :math:`\texttt{INIT_LIQUID_VOLUME}` and :math:`\texttt{CONST_SOLID_VOLUME}`. The bulk liquid volume is a state that evolves with the net flow rate.
2. Constant volume: :math:`\texttt{LIQUID_VOLUME}` and a porosity. The bulk liquid volume is prescribed and does not change.

The porosity is named after what it describes, consistently with the 1D column units.
Particles that have pores hold liquid of their own, so the ratio is the one of the bulk liquid alone and the field is :math:`\texttt{COL_POROSITY}`.
Particles in rapid equilibrium hold no liquid, so all of it is bulk liquid and the field is :math:`\texttt{TOTAL_POROSITY}`.
Without particles, :math:`\texttt{NPARTYPE} = 0`, the field is :math:`\texttt{COL_POROSITY}`.

**The bulk liquid volume is the liquid outside the particles.** The liquid held in the particle pores belongs to the particle volume, not to it.

``INIT_LIQUID_VOLUME``

   Initial volume of the bulk liquid of the variable volume parameterization

   **Unit:** :math:`\mathrm{m}^{3}`

   ================  =========================  =============
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** 1
   ================  =========================  =============

``CONST_SOLID_VOLUME``

   Total volume :math:`V^{p}` of the particles of the variable volume parameterization, **including the liquid in their pores**. Required if :math:`\texttt{NPARTYPE} \geq 1` and the particle types have pores, otherwise optional and defaulting to :math:`0`.

   **Unit:** :math:`\mathrm{m}^{3}`

   ================  =========================  =============
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** 1
   ================  =========================  =============

``LIQUID_VOLUME``

   Prescribed volume of the bulk liquid. Requires the flow rates to cancel in every section, that is :math:`F_{\text{in}} = F_{\text{out}} + F_{\text{filter}}`.

   **Unit:** :math:`\mathrm{m}^{3}`

   ================  ======================  =============
   **Type:** double  **Range:** :math:`> 0`  **Length:** 1
   ================  ======================  =============

``COL_POROSITY`` / ``TOTAL_POROSITY``

   Ratio of the bulk liquid volume to the total volume, :math:`\varepsilon = V^{\ell} / (V^{\ell} + V^{p})`, where :math:`V^{p}` is the total volume of the particles including their pore liquid. Which of the two names applies is described above. Only used together with :math:`\texttt{LIQUID_VOLUME}` and fixes the particle volume to :math:`V^{p} = V^{\ell} (1 - \varepsilon) / \varepsilon`. Required if :math:`\texttt{NPARTYPE} \geq 1`, otherwise optional and defaulting to 1.

   ================  ========================  =============
   **Type:** double  **Range:** :math:`(0,1]`  **Length:** 1
   ================  ========================  =============

Sensitivities are available with respect to :math:`\texttt{CONST_SOLID_VOLUME}` in both parameterizations and with respect to the liquid volume under the name that was used to specify it, that is :math:`\texttt{INIT_LIQUID_VOLUME}` or :math:`\texttt{LIQUID_VOLUME}`.
Note that changing one of the two volumes leaves the other one untouched, so the porosity changes implicitly.

``INIT_STATE``

   Full state vector for initialization (optional, :math:`\texttt{INIT_C}`, :math:`\texttt{INIT_CP}`, :math:`\texttt{INIT_CS}`, and the liquid volume field will be ignored; if length is :math:`2\texttt{NDOF}`, then the second half is used for time derivatives).
   The ordering of the state vector is defined in :ref:`UnitOperationStateOrdering`.

   **Unit:** :math:`various`

   ================  =============================  ====================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NDOF} / 2\texttt{NDOF}`
   ================  =============================  ====================================================

``FLOWRATE_FILTER``

   Flow rate of pure liquid without components (optional, defaults to :math:`0`)

   **Unit:** :math:`\mathrm{m}^{3}\,\mathrm{s}^{-1}`

   ================  =========================  ======================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`1 / \texttt{NSEC}`
   ================  =========================  ======================================

``PAR_TYPE_VOLFRAC``

   Volume fractions of the particle types, have to sum to 1. Required if :math:`\texttt{NPARTYPE} > 1`.

   ================  ========================  =====================================
   **Type:** double  **Range:** :math:`[0,1]`  **Length:** :math:`\texttt{NPARTYPE}`
   ================  ========================  =====================================

``NREAC_LIQUID``

   Number of bulk liquid phase reaction models (optional, only if liquid reactions are present).

   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 0`  **Length:** 1
   =============  =========================  =============


Group /input/model/unit_XXX/particle_type_XXX
----------------------------------------------

One group per particle type, i.e. :math:`\texttt{NPARTYPE}` groups numbered ``particle_type_000``, ``particle_type_001``, and so on.
Each holds the particle geometry, the film and pore transport, the binding model and the particle reactions, see :ref:`particle_model_config`.

A particle type in rapid equilibrium with the bulk liquid only contributes bound states to the state vector, so apart from :math:`\texttt{HAS_FILM_DIFFUSION}`, which is required, the particle transport fields of :ref:`particle_model_config` do not apply to it.
Since the bulk liquid of a well mixed vessel is not flowing, a velocity dependent :math:`\texttt{FILM_DIFFUSION_DEP}` is not supported.

``ADSORPTION_MODEL``

   Specifies the type of binding model of this particle type

   ================  ==========================================  =============
   **Type:** string  **Range:** See Section :ref:`FFAdsorption`  **Length:** 1
   ================  ==========================================  =============

``NBOUND``

   Number of bound states for each component of this particle type. The bound states are kept if :math:`\texttt{ADSORPTION_MODEL}` is :math:`\texttt{NONE}`, which leaves states that solid phase and cross-phase reactions can act on.

   =============  =========================  ==================================
   **Type:** int  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NCOMP}`
   =============  =========================  ==================================

``INIT_CP``

   Initial concentrations for each component in the particle liquid phase of this particle type. This field is optional and defaults to :math:`\texttt{INIT_C}`. A particle type in rapid equilibrium with the bulk liquid has no liquid phase of its own, so specifying it there is an error.

   **Unit:** :math:`\mathrm{mol}\,\mathrm{m}_{\mathrm{MP}}^{-3}`

   ================  =========================  ==================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NCOMP}`
   ================  =========================  ==================================

``INIT_CS``

   Initial concentrations for each bound state of each component in the solid phase of this particle type in component-major ordering. This field is optional and defaults to all 0.

   **Unit:** :math:`\mathrm{mol}\,\mathrm{m}_{\mathrm{SP}}^{-3}`

   ================  =========================  ===================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NBOUND}`
   ================  =========================  ===================================

``NREAC_SOLID``

   Number of solid phase reaction models of this particle type (optional, only if solid phase reactions are present).

   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 0`  **Length:** 1
   =============  =========================  =============

``NREAC_CROSS_PHASE``

   Number of cross-phase reaction models of this particle type (optional, only if cross-phase reactions are present).

   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 0`  **Length:** 1
   =============  =========================  =============


Group /input/model/unit_XXX/particle_type_XXX/adsorption
---------------------------------------------------------
Contains the parameters of the binding model given in :math:`\texttt{ADSORPTION_MODEL}`, see :ref:`FFAdsorption`.


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

The group is optional and additionally holds the settings of the nonlinear solver used for consistent initialization, see :ref:`non_consistency_solver_parameters`.
The spatial discretization of a general rate particle type is given in its own group, see :ref:`particle_model_config`.


Group /input/model/unit_XXX/liquid_reaction_YYY
------------------------------------------------
Each bulk liquid phase reaction is specified in another subgroup `liquid_reaction_YYY`, see :ref:`FFReaction`.
Solid phase and cross-phase reactions are placed in the group of the respective particle type.
