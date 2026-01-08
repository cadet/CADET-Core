.. _particle_model_config:

Particle Model
==============

Group /input/model/unit_XXX/particle_type_XXX
------------------------------------------------------------

``PAR_POROSITY``

   Particle porosity of all particle types or for each particle type
   
   ================  ========================  =============
   **Type:** double  **Range:** :math:`(0,1]`  **Length:** 1
   ================  ========================  =============

``PAR_RADIUS``

   Particle radius of all particle types or for each particle type

   **Unit:** :math:`\mathrm{m}`
   
   ================  =====================  =============
   **Type:** double  **Range:** :math:`>0`  **Length:** 1
   ================  =====================  =============

``PAR_CORERADIUS``

   Particle core radius of all particle types or for each particle type (optional, defaults to :math:`\mathrm{0}`)
   Is only applied when :math:`\texttt{HAS_PORE_DIFFUSION} == 1`.

   **Unit:** :math:`\mathrm{m}`
   
   ================  ===========================================  =============
   **Type:** double  **Range:** :math:`[0, \texttt{PAR_RADIUS})`  **Length:** 1
   ================  ===========================================  =============

``PORE_ACCESSIBILITY``

   Pore accessibility factor of each component in each particle type (optional, defaults to :math:`1`).
   Note: Should not be used in combination with any binding model!
   
   ================  =========================  =============================================================
   **Type:** double  **Range:** :math:`(0, 1]`  **Length:** see :math:`\texttt{PORE_ACCESSIBILITY_MULTIPLEX}`
   ================  =========================  =============================================================

``PORE_ACCESSIBILITY_MULTIPLEX``

   Multiplexing mode of :math:`\texttt{PORE_ACCESSIBILITY}`. Determines whether :math:`\texttt{PORE_ACCESSIBILITY}` is treated as component-, type-, and/or section-independent. This field is optional. When left out, multiplexing behavior is inferred from the length of :math:`\texttt{PORE_ACCESSIBILITY}`. Valid modes are:

   0. Component-dependent, section-independent; length of :math:`\texttt{PORE_ACCESSIBILITY}` is :math:`\texttt{NCOMP}`
   1. Component-dependent, section-dependent; length of :math:`\texttt{PORE_ACCESSIBILITY}` is :math:`\texttt{NCOMP} \cdot \texttt{NSEC}`; ordering is section-major 
   
   =============  ============================  =============
   **Type:** int  **Range:** :math:`\{0, 1 \}`  **Length:** 1
   =============  ============================  =============

``HAS_FILM_DIFFUSION``

	Specifies whether transport into the particles is limited by film diffusion kinetics.

   =============  ============================  =============
   **Type:** int  **Range:** :math:`\{0, 1 \}`  **Length:** 1
   =============  ============================  =============

``FILM_DIFFUSION``

   Film diffusion coefficients for each component of each particle type, required if :math:`\texttt{HAS_FILM_DIFFUSION} == 1`, otherwise ignored.

   **Unit:** :math:`\mathrm{m}\,\mathrm{s}^{-1}`
   
   ================  =========================  =========================================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** see :math:`\texttt{FILM_DIFFUSION_MULTIPLEX}`
   ================  =========================  =========================================================

``FILM_DIFFUSION_MULTIPLEX``

   Multiplexing mode of :math:`\texttt{FILM_DIFFUSION}`. Determines whether :math:`\texttt{FILM_DIFFUSION}` is treated as component-, type-, and/or section-independent.  This field is optional. When left out, multiplexing behavior is inferred from the length of :math:`\texttt{FILM_DIFFUSION}`.  Valid modes are: 

   0. Component-dependent, section-independent; length of :math:`\texttt{FILM_DIFFUSION}` is :math:`\texttt{NCOMP}`
   1. Component-dependent, section-dependent; length of :math:`\texttt{FILM_DIFFUSION}` is :math:`\texttt{NCOMP} \cdot \texttt{NSEC}`; ordering is section-major 
   
   =============  ============================  =============
   **Type:** int  **Range:** :math:`\{0, 1 \}`  **Length:** 1
   =============  ============================  =============

``HAS_PORE_DIFFUSION``

	Specifies whether radial transport within the particle pores is limited by pore diffusion kinetics.

   =============  ============================  =============
   **Type:** int  **Range:** :math:`\{0, 1 \}`  **Length:** 1
   =============  ============================  =============

``PORE_DIFFUSION``

   Effective particle diffusion coefficients of each component in each particle type, required if :math:`\texttt{HAS_PORE_DIFFUSION} == 1`, otherwise ignored.

   **Unit:** :math:`\mathrm{m}_{\mathrm{MP}}^{2}\,\mathrm{s}^{-1}`
   
   ================  ======================  ========================================================
   **Type:** double  **Range:** :math:`> 0`  **Length:** see :math:`\texttt{PORE_DIFFUSION_MULTIPLEX}`
   ================  ======================  ========================================================

``PORE_DIFFUSION_MULTIPLEX``

   Multiplexing mode of :math:`\texttt{PORE_DIFFUSION}`. Determines whether :math:`\texttt{PORE_DIFFUSION}` is treated as component-, type-, and/or section-independent. This field is optional. When left out, multiplexing behavior is inferred from the length of :math:`\texttt{PORE_DIFFUSION}`. Valid modes are: 

   0. Component-dependent, section-independent; length of :math:`\texttt{PORE_DIFFUSION}` is :math:`\texttt{NCOMP}`
   1. Component-dependent, section-dependent; length of :math:`\texttt{PORE_DIFFUSION}` is :math:`\texttt{NCOMP} \cdot \texttt{NSEC}`; ordering is section-major 

   
   =============  ============================  =============
   **Type:** int  **Range:** :math:`\{0, 1 \}`  **Length:** 1
   =============  ============================  =============

``HAS_SURFACE_DIFFUSION``

   Specifies whether radial transport within the particle is supported but limited by surface diffusion kinetics.

   =============  ============================  =============
   **Type:** int  **Range:** :math:`\{0, 1 \}`  **Length:** 1
   =============  ============================  =============

``SURFACE_DIFFUSION``

   Particle surface diffusion coefficients of each bound state of each component in each particle type.
   Ooptional, defaults to all 0 :math:`\mathrm{m}_{\mathrm{SP}}^{2}\,\mathrm{s}^{-1}`, required if :math:`\texttt{HAS_SURFACE_DIFFUSION} == 1`, otherwise ignored.

   **Unit:** :math:`\mathrm{m}_{\mathrm{SP}}^{2}\,\mathrm{s}^{-1}`

   ================  =========================  ============================================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** see :math:`\texttt{SURFACE_DIFFUSION_MULTIPLEX}`
   ================  =========================  ============================================================
   
``SURFACE_DIFFUSION_MULTIPLEX``
   Multiplexing mode of :math:`\texttt{SURFACE_DIFFUSION}`. Determines whether :math:`\texttt{SURFACE_DIFFUSION}` is treated as component-, type-, and/or section-independent. This field is optional. When left out, multiplexing behavior is inferred from the length of :math:`\texttt{SURFACE_DIFFUSION}`. Valid modes are: 

   0. Component-dependent, section-independent; length of :math:`\texttt{SURFACE_DIFFUSION}` is :math:`\texttt{NBOUND}`
   1. Component-dependent, section-dependent; length of :math:`\texttt{SURFACE_DIFFUSION}` is :math:`\texttt{NBOUND} \cdot \texttt{NSEC}`; ordering is section-major 
   
   =============  ============================  =============
   **Type:** int  **Range:** :math:`\{0, 1 \}`  **Length:** 1
   =============  ============================  =============

``PAR_GEOM``

   Specifies the particle geometry for all or each particle type.
   Valid values are :math:`\texttt{SPHERE}`, :math:`\texttt{CYLINDER}`, :math:`\texttt{SLAB}`. Optional, defaults to :math:`\texttt{SPHERE}`.
   Is only applied when :math:`\texttt{HAS_PORE_DIFFUSION} == 1`.
   
   ================  =========================================================================  =============
   **Type:** string  **Range:** :math:`\{\texttt{SPHERE}, \texttt{CYLINDER}, \texttt{SLAB} \}`  **Length:** 1
   ================  =========================================================================  =============


``ADSORPTION_MODEL``

   Specifies the type of binding model of each particle type (or of all particle types if length is :math:`1`)
   
   ================  ==============================  =============
   **Type:** string  **Range:** :ref:`FFAdsorption`  **Length:** 1
   ================  ==============================  =============

``NBOUND``

   Number of bound states for each component in each particle type in particle type major ordering
   
   =============  =========================  ==========================================
   **Type:** int  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NCOMP}`
   =============  =========================  ==========================================

``REACTION_MODEL``

   Specifies the type of reaction model of each particle type (or of all particle types if length is :math:`1`).
   The model is configured in the subgroup :math:`\texttt{reaction_particle}`, or :math:`\texttt{reaction_particle_XXX}` in case of disabled multiplexing.
   
   ================  ========================================  =============
   **Type:** string  **Range:** See Section :ref:`FFReaction`  **Length:** 1
   ================  ========================================  =============

``INIT_CP``

   Initial concentrations for each component in the bead liquid phase (optional, :math:`\texttt{INIT_C}` is used if left out).
   The length of this field is :math:`\texttt{NCOMP}`.
   Only in case of a 2D bulk model, the field length *can* also be :math:`\texttt{NCOMP} \cdot \texttt{NRAD}` to specify radial dependence, in which case the ordering is radial position major.
   If :math:`\texttt{BINDING_PARTYPE_DEPENDENT}` is :math:`0`, the values across different particle types (if multiple particle types are being used) must be the same.

   **Unit:** :math:`\mathrm{mol}\,\mathrm{m}_{\mathrm{MP}}^{-3}`
   
   ================  =========================  ==================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NCOMP}`
   ================  =========================  ==================================

``INIT_CS``

   Initial concentrations for each bound state of each component in the bead solid phase.
   The length of this field is :math:`\texttt{NBOUND}`.
   Only in case of a 2D bulk model, the field length *can* also be :math:`\texttt{NBOUND} \cdot \texttt{NRAD}` to specify radial dependence, in which case the ordering is radial position major.
   If :math:`\texttt{BINDING_PARTYPE_DEPENDENT}` is :math:`0`, the values across different particle types (if multiple particle types are being used) must be the same.

   **Unit:** :math:`\mathrm{mol}\,\mathrm{m}_{\mathrm{SP}}^{-3}`
   
   ================  =========================
   **Type:** double  **Range:** :math:`\geq 0`
   ================  =========================

``PARAMNAME_PARTYPE_DEPENDENT``
	
	Only required for parameter sensitivities of the respective parameter, defaults to 1.
	Specifies whether or not a parameter (parameter name substitutes the `PARAMNAME` prefix of field name) is the same across particle types or 'dependent' on the particle type.
	In the first case, the parameter sensitivity should be computed jointly across particle types, by specifying this field as 0.
	Can be specified for any of the above parameters except `NCOMP` and `NBOUND`.
	For more information on parameter sensitivities, see :ref:`spatial_discretization_methods`.
   
   =============  ============================  =============
   **Type:** int  **Range:** :math:`\{0, 1 \}`  **Length:** 1
   =============  ============================  =============

Group /input/model/unit_XXX/particle_type_XXX/discretization
------------------------------------------------------------

``PAR_DISC_TYPE``

   Specifies the discretization scheme inside the particles for all or each particle type. Valid values are :math:`\texttt{EQUIDISTANT}`, :math:`\texttt{EQUIVOLUME}`, and :math:`\texttt{USER_DEFINED}`.
   
   ================  =================================================
   **Type:** string  **Length:** 1
   ================  =================================================

``PAR_DISC_VECTOR``

   Node coordinates for the DG element (or FV cell) boundaries (ignored if :math:`\texttt{PAR_DISC_TYPE} \neq \texttt{USER_DEFINED}`).
   The coordinates are relative and have to include the endpoints :math:`0` and :math:`1`.
   They are later linearly mapped to the true radial range :math:`[r_{c,j}, r_{p,j}]`.
   The coordinates for each particle type are appended to one long vector in particle type major ordering.
   
   ================  ========================  ======================================
   **Type:** double  **Range:** :math:`[0,1]`  **Length:** :math:`\texttt{NELEM} + 1`
   ================  ========================  ======================================

Spatial discretization - Numerical Methods
------------------------------------------

CADET offers two spatial discretization methods: Finite Volumes (FV) and Discontinuous Galerkin (DG). Each method has it's own set of input fields.
While both methods approximate the same solution to the same underlying model, they may differ in terms of computational performance.
With our currently implemented variants of FV and DG, FV perform better for solutions with steep gradients or discontinuities, while DG can be much faster for rather smooth solutions.
For the same number of discrete points, DG will generally be slower but often more accurate.

For further information on the choice of discretization methods and their parameters, see :ref:`sensitivity`.

``SPATIAL_METHOD``

   Spatial discretization method. Optional, defaults to :math:`\texttt{DG}`

   ================  ===============================================  =============
   **Type:** string  **Range:** :math:`\{\texttt{FV}, \texttt{DG}\}`  **Length:** 1
   ================  ===============================================  =============

Discontinuous Galerkin
----------------------

``PAR_POLYDEG``

   DG particle (radial) polynomial degree. Optional, defaults to 3. The total number of particle (radial) discrete points is given by (``PARPOLYDEG`` + 1 ) * ``PAR_NELEM``.
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``PAR_NELEM``

   Number of particle (radial) discretization DG elements for each particle type. For the particle discretization, it is usually most performant to fix ``PAR_NELEM`` = 1 and to increase the polynomial degree for more accuracy.
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``PAR_GSM``

   Specifies whether Galerkin spectral method should be used (as opposed to discontinuous variant, DGSEM), optional, defaults to 1 if ``PAR_NELEM`` == 1. Always recommended.
   
   =============  ============================  =============
   **Type:** int  **Range:** :math:`\{ 0,1 \}`  **Length:** 1
   =============  ============================  =============
   
   For further discretization parameters, see also :ref:`non_consistency_solver_parameters`.

Finite Volumes
--------------

``NCELLS``

   Number of particle (radial) discretization points for each particle type
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``FV_BOUNDARY_ORDER``

   Order of accuracy of outer particle boundary condition. Optional, defaults to :math:`2`.
   
   =============  ============================  =============
   **Type:** int  **Range:** :math:`\{ 1,2 \}`  **Length:** 1
   =============  ============================  =============

``OPTIMIZE_PAR_BANDWIDTH``

   Only available if both particle and bulk spatial methods are all FV
   Determines whether the surface diffusion parameters :math:`\texttt{SURFACE_DIFFUSION}` are fixed if the parameters are zero.
   If the parameters are fixed to zero (:math:`\texttt{FIX_ZERO_SURFACE_DIFFUSION} = 1`, :math:`\texttt{SURFACE_DIFFUSION} = 0`), the parameters must not become non-zero during this or subsequent simulation runs.
   The internal data structures are optimized for a more efficient simulation.  This field is optional and defaults to :math:`0` (optimization disabled in favor of flexibility).
   
   =============  ===========================  =============
   **Type:** int  **Range:** :math:`\{0, 1\}`  **Length:** 1
   =============  ===========================  =============

When using the FV method, we generally recommend specifying ``USE_MODIFIED_NEWTON = 0`` in :ref:`FFSolverTime`, i.e. to use the full Newton method to solve the linear system within the time integrator.
For further discretization parameters, see also :ref:`flux_reconstruction_methods` (FV specific)), and :ref:`non_consistency_solver_parameters`.
