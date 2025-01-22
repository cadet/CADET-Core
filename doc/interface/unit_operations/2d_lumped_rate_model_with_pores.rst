.. _2d_lumped_rate_model_with_pores_config:

Two dimensional lumped rate model with pores
============================================

Group /input/model/unit_XXX - UNIT_TYPE = LUMPED_RATE_MODEL_WITH_PORES_2D
-------------------------------------------------------------------------

For information on model equations, refer to :ref:`2d_lumped_rate_model_with_pores_model`.


``UNIT_TYPE``

   Specifies the type of unit operation model
   
   ================  ===========================================================  =============
   **Type:** string  **Range:** :math:`\texttt{LUMPED_RATE_MODEL_WITH_PORES_2D}`  **Length:** 1
   ================  ===========================================================  =============
   
``NCOMP``

   Number of chemical components in the chromatographic medium
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============
   
``ADSORPTION_MODEL``

   Specifies the type of binding model of each particle type (or of all particle types if length is 1)
   
   ================  ==========================================  ==========================================
   **Type:** string  **Range:** See Section :ref:`FFAdsorption`  **Length:** :math:`1 / \texttt{NPARTYPE}`
   ================  ==========================================  ==========================================
   
``ADSORPTION_MODEL_MULTIPLEX``

   Multiplexing mode of :math:`\texttt{ADSORPTION_MODEL}`. If set to 0, each particle type has a different binding model and the length of :math:`\texttt{ADSORPTION_MODEL}` is :math:`\texttt{NPARTYPE}`. If set to 1, all particle types share the same binding model and the length of :math:`\texttt{ADSORPTION_MODEL}` is 1.  This field is optional and inferred from the length of :math:`\texttt{ADSORPTION_MODEL}` if left out.
   
   =============  ===========================  =============
   **Type:** int  **Range:** :math:`\{0, 1\}`  **Length:** 1
   =============  ===========================  =============
   
``NBOUND``

   Number of bound states for each component in each particle type in type-major ordering
   
   =============  =========================  ==========================================================================
   **Type:** int  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NCOMP} / \texttt{NPARTYPE} \cdot \texttt{NCOMP}`
   =============  =========================  ==========================================================================
   
``REACTION_MODEL_BULK``

   Specifies the type of reaction model of the bulk volume. The model is configured in the subgroup :math:`\texttt{reaction_bulk}`.
   
   ================  ========================================  =============
   **Type:** string  **Range:** See Section :ref:`FFReaction`  **Length:** 1
   ================  ========================================  =============
   
``REACTION_MODEL_PARTICLES``

   Specifies the type of reaction model of each particle type (or of all particle types if length is 1). The model is configured in the subgroup :math:`\texttt{reaction_particle}`, or :math:`\texttt{reaction_particle_XXX}` in case of disabled multiplexing.
   
   ================  ========================================  ==========================================
   **Type:** string  **Range:** See Section :ref:`FFReaction`  **Length:** :math:`1 / \texttt{NPARTYPE}`
   ================  ========================================  ==========================================
   
``REACTION_MODEL_PARTICLES_MULTIPLEX``

   Multiplexing mode of :math:`\texttt{REACTION_MODEL_PARTICLES}`. If set to 0, each particle type has a different reaction model and the length of :math:`\texttt{REACTION_MODEL_PARTICLES}` is :math:`\texttt{NPARTYPE}`. If set to 1, all particle types share the same reaction model and the length of :math:`\texttt{REACTION_MODEL_PARTICLES}` is 1.  This field is optional and inferred from the length of :math:`\texttt{REACTION_MODEL_PARTICLES}` if left out.
   
   =============  ===========================  =============
   **Type:** int  **Range:** :math:`\{0, 1\}`  **Length:** 1
   =============  ===========================  =============
   
``INIT_C``

   Initial concentrations for each component in all radial zones the bulk mobile phase (length :math:`\texttt{NCOMP}`), or for each component in each radial zone (length :math:`\texttt{NCOMP} \cdot \texttt{RAD_NODES}`, ordering radial-major)

   **Unit:** :math:`\mathrm{mol}\,\mathrm{m}_{\mathrm{IV}}^{-3}`
   
   ================  =========================  ============================================================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NCOMP} / \texttt{NCOMP} \cdot \texttt{RAD_NODES}`
   ================  =========================  ============================================================================
   
``INIT_CP``

   Initial concentrations for each component in the bead liquid phase (optional, :math:`\texttt{INIT_C}` is used if left out). The length of this field can be :math:`\texttt{NCOMP}` (same values for each radial zone and particle type), :math:`\texttt{NPARTYPE} \cdot \texttt{NCOMP}` (same values for each radial zone), :math:`\texttt{RAD} \cdot \texttt{NCOMP}` (same values for each particle type), or :math:`\texttt{RAD_NODES} \cdot \texttt{NPARTYPE} \cdot \texttt{NCOMP}`. The ordering is radial-type-major.  Values for each particle type can only be given when :math:`\texttt{ADSORPTION_MODEL_MULTIPLEX}` is 0. In the radial-inhomogeneous case, the :math:`\texttt{SENS_REACTION}` field is used for indexing the radial zone when specifying parameter sensitivities.

   **Unit:** :math:`\mathrm{mol}\,\mathrm{m}_{\mathrm{MP}}^{-3}`
   
   ================  =========================
   **Type:** double  **Range:** :math:`\geq 0`
   ================  =========================
   
``INIT_Q``

   Initial concentrations for each bound state of each component in the bead solid phase. If :math:`\texttt{ADSORPTION_MODEL_MULTIPLEX}` is 0, values for each particle type are required in type-component-major ordering (length is :math:`\texttt{NTOTALBND}`). If :math:`\texttt{ADSORPTION_MODEL_MULTIPLEX}` is 1, values for one particle type are required in component-major ordering (length is :math:`\sum_{i = 0}^{\texttt{NCOMP} - 1} \texttt{NBND}_i`).  Alternatively, values for each radial zone can be supplied. If :math:`\texttt{ADSORPTION_MODEL_MULTIPLEX}` is 0, values for each radial zone and each particle type are required in radial-type-component-major ordering (length is :math:`\texttt{RAD_NODES} \cdot \texttt{NTOTALBND}`). If :math:`\texttt{ADSORPTION_MODEL_MULTIPLEX}` is 1, values for each radial zone and all particle types are required in radial-component-major ordering (length is :math:`\texttt{RAD_NODES} \cdot \sum_{i = 0}^{\texttt{NCOMP} - 1} \texttt{NBND}_i`). In the radial-inhomogeneous case, the :math:`\texttt{SENS_REACTION}` field is used for indexing the radial zone when specifying parameter sensitivities.

   **Unit:** :math:`\mathrm{mol}\,\mathrm{m}_{\mathrm{SP}}^{-3}`
   
   ================  =========================
   **Type:** double  **Range:** :math:`\geq 0`
   ================  =========================
   
``INIT_STATE``

   Full state vector for initialization (optional, :math:`\texttt{INIT_C}`, :math:`\texttt{INIT_CP}`, and :math:`\texttt{INIT_Q}` will be ignored; if length is :math:`2\texttt{NDOF}`, then the second half is used for time derivatives)

   **Unit:** :math:`various`
   
   ================  =============================  ==================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NDOF} / 2\texttt{NDOF}`
   ================  =============================  ==================================================
   
``COL_DISPERSION``

   Axial dispersion coefficient.  In case of a spatially inhomogeneous setting, the :math:`\texttt{SENS_PARTYPE}` field is used for indexing the radial zone when specifying parameter sensitivities.

   **Unit:** :math:`\mathrm{m}_{\mathrm{IV}}^{2}\,\mathrm{s}^{-1}`
   
   ================  =========================  =========================================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** see :math:`\texttt{COL_DISPERSION_MULTIPLEX}`
   ================  =========================  =========================================================
   
``COL_DISPERSION_MULTIPLEX``

   Multiplexing mode of :math:`\texttt{COL_DISPERSION}`. Determines whether :math:`\texttt{COL_DISPERSION}` is treated as component-, radial-, and/or section-independent.  This field is optional. When left out, multiplexing behavior is inferred from the length of :math:`\texttt{COL_DISPERSION}`.  Valid modes are: 

  0. Component-independent, radial-independent, section-independent; length of :math:`\texttt{COL_DISPERSION}` is 1 
  1. Component-independent, radial-dependent, section-independent; length of :math:`\texttt{COL_DISPERSION}` is :math:`\texttt{RAD_NELEM}` 
  2. Component-dependent, radial-independent, section-independent; length of :math:`\texttt{COL_DISPERSION}` is :math:`\texttt{NCOMP}` 
  3. Component-dependent, radial-dependent, section-independent; length of :math:`\texttt{COL_DISPERSION}` is :math:`\texttt{NCOMP} \cdot \texttt{RAD_NELEM}`; ordering is radial-major 
  4. Component-independent, radial-independent, section-dependent; length of :math:`\texttt{COL_DISPERSION}` is :math:`\texttt{NSEC}` 
  5. Component-independent, radial-dependent, section-dependent; length of :math:`\texttt{COL_DISPERSION}` is :math:`\texttt{RAD_NELEM} \cdot \texttt{NSEC}`; ordering is section-major 
  6. Component-dependent, radial-independent, section-independent; length of :math:`\texttt{COL_DISPERSION}` is :math:`\texttt{NCOMP} \cdot \texttt{NSEC}`; ordering is section-major 
  7. Component-dependent, radial-dependent, section-dependent; length of :math:`\texttt{COL_DISPERSION}` is :math:`\texttt{NCOMP} \cdot \texttt{RAD_NELEM} \cdot \texttt{NSEC}`; ordering is section-radial-major 
   
   =============  ===================================  =============
   **Type:** int  **Range:** :math:`\{0, \dots, 7 \}`  **Length:** 1
   =============  ===================================  =============
   
``COL_DISPERSION_RADIAL``

   Radial dispersion coefficient.  In case of a spatially inhomogeneous setting, the :math:`\texttt{SENS_PARTYPE}` field is used for indexing the radial zone when specifying parameter sensitivities.

   **Unit:** :math:`\mathrm{m}_{\mathrm{IV}}^{2}\,\mathrm{s}^{-1}`
   
   ================  =========================  ================================================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** see :math:`\texttt{COL_DISPERSION_RADIAL_MULTIPLEX}`
   ================  =========================  ================================================================
   
``COL_DISPERSION_RADIAL_MULTIPLEX``

   Multiplexing mode of :math:`\texttt{COL_DISPERSION_RADIAL}`. Determines whether :math:`\texttt{COL_DISPERSION_RADIAL}` is treated as component-, radial-, and/or section-independent.  This field is optional. When left out, multiplexing behavior is inferred from the length of :math:`\texttt{COL_DISPERSION_RADIAL}`.  Valid modes are: 

  0. Component-independent, radial-independent, section-independent; length of :math:`\texttt{COL_DISPERSION_RADIAL}` is 1 
  1. Component-independent, radial-dependent, section-independent; length of :math:`\texttt{COL_DISPERSION_RADIAL}` is :math:`\texttt{RAD_NELEM}` 
  2. Component-dependent, radial-independent, section-independent; length of :math:`\texttt{COL_DISPERSION_RADIAL}` is :math:`\texttt{NCOMP}` 
  3. Component-dependent, radial-dependent, section-independent; length of :math:`\texttt{COL_DISPERSION_RADIAL}` is :math:`\texttt{NCOMP} \cdot \texttt{RAD_NELEM}`; ordering is radial-major 
  4. Component-independent, radial-independent, section-dependent; length of :math:`\texttt{COL_DISPERSION_RADIAL}` is :math:`\texttt{NSEC}` 
  5. Component-independent, radial-dependent, section-dependent; length of :math:`\texttt{COL_DISPERSION_RADIAL}` is :math:`\texttt{RAD_NELEM} \cdot \texttt{NSEC}`; ordering is section-major 
  6. Component-dependent, radial-independent, section-independent; length of :math:`\texttt{COL_DISPERSION_RADIAL}` is :math:`\texttt{NCOMP} \cdot \texttt{NSEC}`; ordering is section-major 
  7. Component-dependent, radial-dependent, section-dependent; length of :math:`\texttt{COL_DISPERSION_RADIAL}` is :math:`\texttt{NCOMP} \cdot \texttt{RAD_NELEM} \cdot \texttt{NSEC}`; ordering is section-radial-major 
   
   =============  ===================================  =============
   **Type:** int  **Range:** :math:`\{0, \dots, 7 \}`  **Length:** 1
   =============  ===================================  =============
   
``COL_LENGTH``

   Column length

   **Unit:** :math:`\mathrm{m}`
   
   ================  ======================  =============
   **Type:** double  **Range:** :math:`> 0`  **Length:** 1
   ================  ======================  =============
   
``COL_RADIUS``

   Column radius. This parameter is optional if ``CROSS_SECTION_AREA`` is provided.

   **Unit:** :math:`\mathrm{m}`

   ================  ======================  =============
   **Type:** double  **Range:** :math:`> 0`  **Length:** 1
   ================  ======================  =============

``CROSS_SECTION_AREA``

   Cross section area of the column. This parameter is optional and will be ignored if `COL_RADIUS` is provided.

   **Unit:** :math:`\mathrm{m}^{2}`

   ================  =====================  =============
   **Type:** double  **Range:** :math:`>0`  **Length:** 1
   ================  =====================  =============
   
``COL_POROSITY``

   Column porosity, either constant (length is 1) or for each radial zone (length is :math:`\texttt{RAD_NELEM}`).  In case of a spatially inhomogeneous setting, the :math:`\texttt{SENS_PARTYPE}` field is used for indexing the radial zone when specifying parameter sensitivities.
   
   ================  ========================  =====================================
   **Type:** double  **Range:** :math:`(0,1]`  **Length:** :math:`1 / \texttt{RAD_NELEM}`
   ================  ========================  =====================================
   
``FILM_DIFFUSION``

   Film diffusion coefficients for each component of each particle type

   **Unit:** :math:`\mathrm{m}\,\mathrm{s}^{-1}`
   
   ================  =========================  =======================================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** see :math:`\texttt{FILM_DIFFUSION_MULTIPLEX}`
   ================  =========================  =======================================================
   
``FILM_DIFFUSION_MULTIPLEX``

   Multiplexing mode of :math:`\texttt{FILM_DIFFUSION}`. Determines whether :math:`\texttt{FILM_DIFFUSION}` is treated as component-, type-, and/or section-independent.  This field is optional. When left out, multiplexing behavior is inferred from the length of :math:`\texttt{FILM_DIFFUSION}`.  Valid modes are: 

  0. Component-dependent, type-independent, section-independent; length of :math:`\texttt{FILM_DIFFUSION}` is :math:`\texttt{NCOMP}` 
  1. Component-dependent, type-independent, section-dependent; length of :math:`\texttt{FILM_DIFFUSION}` is :math:`\texttt{NCOMP} \cdot \texttt{NSEC}`; ordering is section-major 
  2. Component-dependent, type-dependent, section-independent; length of :math:`\texttt{FILM_DIFFUSION}` is :math:`\texttt{NCOMP} \cdot \texttt{NPARTYPE}`; ordering is type-major 
  3. Component-dependent, type-dependent, section-dependent; length of :math:`\texttt{FILM_DIFFUSION}` is :math:`\texttt{NCOMP} \cdot \texttt{NPARTYPE} \cdot \texttt{NSEC}`; ordering is section-type-major 
   
   =============  ===================================  =============
   **Type:** int  **Range:** :math:`\{0, \dots, 3 \}`  **Length:** 1
   =============  ===================================  =============
   
``PAR_POROSITY``

   Particle porosity of all particle types or for each particle type
   
   ================  ========================  =========================================
   **Type:** double  **Range:** :math:`(0,1]`  **Length:** :math:`1 / \texttt{NPARTYPE}`
   ================  ========================  =========================================
   
``PAR_RADIUS``

   Particle radius of all particle types or for each particle type

   **Unit:** :math:`\mathrm{m}`
   
   ================  =====================  =========================================
   **Type:** double  **Range:** :math:`>0`  **Length:** :math:`1 / \texttt{NPARTYPE}`
   ================  =====================  =========================================
   
``PORE_ACCESSIBILITY``

   Pore accessibility factor of each component in each particle type (optional, defaults to 1).
   Note: Should not be used in combination with any binding model!
   
   ================  =========================  =============================================================
   **Type:** double  **Range:** :math:`(0, 1]`  **Length:** see :math:`\texttt{PORE_ACCESSIBILITY_MULTIPLEX}`
   ================  =========================  =============================================================
   
``PORE_ACCESSIBILITY_MULTIPLEX``

   Multiplexing mode of :math:`\texttt{PORE_ACCESSIBILITY}`. Determines whether :math:`\texttt{PORE_ACCESSIBILITY}` is treated as component-, type-, and/or section-independent.  This field is optional. When left out, multiplexing behavior is inferred from the length of :math:`\texttt{PORE_ACCESSIBILITY}`.  Valid modes are: 

  0. Component-dependent, type-independent, section-independent; length of :math:`\texttt{PORE_ACCESSIBILITY}` is :math:`\texttt{NCOMP}` 
  1. Component-dependent, type-independent, section-dependent; length of :math:`\texttt{PORE_ACCESSIBILITY}` is :math:`\texttt{NCOMP} \cdot \texttt{NSEC}`; ordering is section-major 
  2. Component-dependent, type-dependent, section-independent; length of :math:`\texttt{PORE_ACCESSIBILITY}` is :math:`\texttt{NCOMP} \cdot \texttt{NPARTYPE}`; ordering is type-major 
  3. Component-dependent, type-dependent, section-dependent; length of :math:`\texttt{PORE_ACCESSIBILITY}` is :math:`\texttt{NCOMP} \cdot \texttt{NPARTYPE} \cdot \texttt{NSEC}`; ordering is section-type-major 
   
   =============  ===================================  =============
   **Type:** int  **Range:** :math:`\{0, \dots, 3 \}`  **Length:** 1
   =============  ===================================  =============
   
``VELOCITY``

   Indicates flow direction in each radial zone (forward if value is positive, backward if value is negative), see Section :ref:`MUOPGRMflow2D`).  In case of a spatially inhomogeneous setting, the :math:`\texttt{SENS_PARTYPE}` field is used for indexing the radial cell when specifying parameter sensitivities.
   
   ================  =============================  ===================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** see :math:`\texttt{VELOCITY_MULTIPLEX}`
   ================  =============================  ===================================================
   
``VELOCITY_MULTIPLEX``

   Multiplexing mode of :math:`\texttt{VELOCITY}`. Determines whether :math:`\texttt{VELOCITY}` is treated as radial- and/or section-independent.  This field is optional. When left out, multiplexing behavior is inferred from the length of :math:`\texttt{VELOCITY}`.  Valid modes are: 

  0. Radial-independent, section-independent; length of :math:`\texttt{VELOCITY}` is 1 
  1. Radial-dependent, section-independent; length of :math:`\texttt{VELOCITY}` is :math:`\texttt{RAD_NELEM}` 
  2. Section-dependent; length of :math:`\texttt{VELOCITY}` is :math:`\texttt{NSEC}` 
  3. Radial-dependent, section-dependent; length of :math:`\texttt{VELOCITY}` is :math:`\texttt{RAD_NELEM} \cdot \texttt{NSEC}`; ordering is section-major 
   
   =============  ===================================  =============
   **Type:** int  **Range:** :math:`\{0, \dots, 3 \}`  **Length:** 1
   =============  ===================================  =============
   
``NPARTYPE``

   Number of particle types. Optional, inferred from the length of :math:`\texttt{NPAR}` or :math:`\texttt{NBOUND}` if left out.
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============
   
``PAR_TYPE_VOLFRAC``

   Volume fractions of the particle types. The volume fractions can be set homogeneous or individually along both axes. For each cell, the volume fractions have to sum to 1.  In case of a spatially inhomogeneous setting, the :math:`\texttt{SENS_SECTION}` field is used for indexing the axial cell and the :math:`\texttt{SENS_REACTION}` field is used for indexing the radial cell when specifying parameter sensitivities.  This field is optional in case of only one particle type.
   
   ================  ========================  ===========================================================
   **Type:** double  **Range:** :math:`[0,1]`  **Length:** see :math:`\texttt{PAR_TYPE_VOLFRAC_MULTIPLEX}`
   ================  ========================  ===========================================================
   
``PAR_TYPE_VOLFRAC_MULTIPLEX``

   Multiplexing mode of :math:`\texttt{PAR_TYPE_VOLFRAC}`. Determines whether :math:`\texttt{PAR_TYPE_VOLFRAC}` is treated as radial- and/or section-independent.  This field is optional. When left out, multiplexing behavior is inferred from the length of :math:`\texttt{PAR_TYPE_VOLFRAC}`.  Valid modes are: 

  0. Radial-independent, axial-independent; length of :math:`\texttt{PAR_TYPE_VOLFRAC}` is :math:`\texttt{NPARTYPE}` 
  1. Radial-dependent, axial-independent; length of :math:`\texttt{PAR_TYPE_VOLFRAC}` is :math:`\texttt{RAD_NELEM} \cdot \texttt{NPARTYPE}`; ordering is radial-major 
  2. Axial-dependent; length of :math:`\texttt{PAR_TYPE_VOLFRAC}` is :math:`\texttt{NCOL} \cdot \texttt{NPARTYPE}`; ordering is axial-major 
  3. Radial-dependent, axial-dependent; length of :math:`\texttt{PAR_TYPE_VOLFRAC}` is :math:`\texttt{NCOL} \cdot \texttt{RAD_NELEM} \cdot \texttt{NPARTYPE}`; ordering is axial-radial-major 
   
   =============  ===================================  =============
   **Type:** int  **Range:** :math:`\{0, \dots, 3 \}`  **Length:** 1
   =============  ===================================  =============
   

Group /input/model/unit_XXX/discretization - UNIT_TYPE - LUMPED_RATE_MODEL_WITH_PORES_2D
----------------------------------------------------------------------------------------

For the spatial discretization of the two-dimensional LRMP, CADET utilizes a Discontinuous Galerkin (DG) method.

For further information on the choice of discretization methods and their parameters, see :ref:`spatial_discretization_methods

``AX_POLYDEG``

   DG polynomial degree in axial direction. Optional, defaults to 4 and :math:`N_d \in \{3, 4, 5\}` is recommended. The total number of axial discrete points is given by (``AX_POLYDEG`` + 1 ) * ``AX_NELEM``
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``AX_NELEM``

   Number of axial DG cells\elements. The total number of axial discrete points is given by (``AX_POLYDEG`` + 1 ) * ``AX_NELEM``
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============
   
``RAD_POLYDEG``

   DG polynomial degree in radial direction. Optional, defaults to 4 and :math:`N_d \in \{3, 4, 5\}` is recommended. The total number of radial discrete points is given by (``RAD_POLYDEG`` + 1 ) * ``RAD_NELEM``
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``RAD_NELEM``

   Number of radial DG cells\elements. The total number of radial discrete points is given by ``RAD_NODES`` = (``RAD_POLYDEG`` + 1 ) * ``RAD_NELEM``
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============
   
``RADIAL_DISC_TYPE``

   Specifies the radial discretization scheme. Valid values are :math:`\texttt{EQUIDISTANT}`, :math:`\texttt{EQUIVOLUME}`, and :math:`\texttt{USER_DEFINED}`.
   
   ================  =============
   **Type:** string  **Length:** 1
   ================  =============
   
``RADIAL_COMPARTMENTS``

   Coordinates for the radial compartment boundaries (ignored if :math:`\texttt{RADIAL_DISC_TYPE} \neq \texttt{USER_DEFINED}`). The coordinates are absolute and have to include the endpoints 0 and :math:`\texttt{COLUMN_RADIUS}`. The values are expected in ascending order (i.e., 0 is the first and :math:`\texttt{COLUMN_RADIUS}` the last value in the array).

   **Unit:** :math:`\mathrm{m}`
   
   ================  =============================================  ==========================================
   **Type:** double  **Range:** :math:`[0,\texttt{COLUMN_RADIUS}]`  **Length:** :math:`\texttt{RAD_NELEM} + 1`
   ================  =============================================  ==========================================
   
   
``USE_ANALYTIC_JACOBIAN``

   Determines whether analytically computed Jacobian matrix (faster) is used (value is 1) instead of Jacobians generated by algorithmic differentiation (slower, value is 0)
   
   =============  ===========================  =============
   **Type:** int  **Range:** :math:`\{0, 1\}`  **Length:** 1
   =============  ===========================  =============
   
