.. _axial_flow_column_1D_config:

Axial Flow Column 1D
====================

Group /input/model/unit_XXX - UNIT_TYPE - COLUMN_MODEL_1D
------------------------------------------------------------

``UNIT_TYPE``

   Specifies the type of unit operation model
   
   ================  ===========================================  =============
   **Type:** string  **Range:** :math:`\texttt{COLUMN_MODEL_1D}`  **Length:** 1
   ================  ===========================================  =============

``NCOMP``

   Number of chemical components in the chromatographic medium
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``CROSS_SECTION_AREA``

   Cross section area of the column (optional if :math:`\texttt{VELOCITY}` is present, see Section :ref:`MUOPGRMflow`)
   **Unit:** :math:`\mathrm{m}^{2}`
   
   ================  =====================  =============
   **Type:** double  **Range:** :math:`>0`  **Length:** 1
   ================  =====================  =============

``COL_LENGTH``

   Column length

   **Unit:** :math:`\mathrm{m}`
   
   ================  ======================  =============
   **Type:** double  **Range:** :math:`> 0`  **Length:** 1
   ================  ======================  =============

``COL_POROSITY``

   Column porosity
   
   ================  ========================  =============
   **Type:** double  **Range:** :math:`(0,1]`  **Length:** 1
   ================  ========================  =============

``NPARTYPE``

   Number of particle types.
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``PAR_TYPE_VOLFRAC``

   Volume fractions of the particle types. The volume fractions can be set for all axial cells together or for each individual axial cell. For each cell, the volume fractions have to sum to :math:`1`. In case of a spatially inhomogeneous setting, the data is expected in cell-major ordering and the :math:`\texttt{SENS_SECTION}` field is used for indexing the axial cell when specifying parameter sensitivities.  This field is optional in case of only one particle type.
   
   ================  ========================  =============================================================================
   **Type:** double  **Range:** :math:`[0,1]`  **Length:** :math:`\texttt{NPARTYPE} / \texttt{NCOL} \cdot \texttt{NPARTYPE}`
   ================  ========================  =============================================================================

``VELOCITY``

   Interstitial velocity of the mobile phase (optional if :math:`\texttt{CROSS_SECTION_AREA}` is present, see Section :ref:`MUOPGRMflow`)
   **Unit:** :math:`\mathrm{m}\,\mathrm{s}^{-1}`
   
   ================  =============================  =======================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`1 / \texttt{NSEC}`
   ================  =============================  =======================================

``COL_DISPERSION``

   Axial dispersion coefficient

   **Unit:** :math:`\mathrm{m}_{\mathrm{IV}}^{2}\,\mathrm{s}^{-1}`
   
   ================  =========================  =========================================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** see :math:`\texttt{COL_DISPERSION_MULTIPLEX}`
   ================  =========================  =========================================================

``COL_DISPERSION_MULTIPLEX``

   Multiplexing mode of :math:`\texttt{COL_DISPERSION}`. Determines whether :math:`\texttt{COL_DISPERSION}` is treated as component- and/or section-independent.  This field is optional. When left out, multiplexing behavior is inferred from the length of :math:`\texttt{COL_DISPERSION}`.  Valid modes are: 

   0. Component-independent, section-independent; length of :math:`\texttt{COL_DISPERSION}` is :math:`1` 
   1. Component-dependent, section-independent; length of :math:`\texttt{COL_DISPERSION}` is :math:`\texttt{NCOMP}` 
   2. Component-independent, section-dependent; length of :math:`\texttt{COL_DISPERSION}` is :math:`\texttt{NSEC}` 
   3. Component-dependent, section-dependent; length of :math:`\texttt{COL_DISPERSION}` is :math:`\texttt{NCOMP} \cdot \texttt{NSEC}`; ordering is section-major 
   
   =============  ===================================  =============
   **Type:** int  **Range:** :math:`\{0, \dots, 3 \}`  **Length:** 1
   =============  ===================================  =============

``REACTION_MODEL_BULK``

   Specifies the type of reaction model of the bulk volume. The model is configured in the subgroup :math:`\texttt{reaction_bulk}`.
   
   ================  ========================================  =============
   **Type:** string  **Range:** See Section :ref:`FFReaction`  **Length:** 1
   ================  ========================================  =============

``INIT_C``

   Initial concentrations for each component in the bulk mobile phase

   **Unit:** :math:`\mathrm{mol}\,\mathrm{m}_{\mathrm{IV}}^{-3}`
   
   ================  =========================  ==================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NCOMP}`
   ================  =========================  ==================================

``INIT_STATE``

   Full state vector for initialization (optional, :math:`\texttt{INIT_C}`, :math:`\texttt{INIT_CP}`, and :math:`\texttt{INIT_CS}` will be ignored; if length is :math:`2\texttt{NDOF}`, then the second half is used for time derivatives).
   The ordering of the state vector is defined in :ref:`UnitOperationStateOrdering`.

   **Unit:** :math:`various`
   
   ================  =============================  ==================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NDOF} / 2\texttt{NDOF}`
   ================  =============================  ==================================================


Group /input/model/unit_XXX/particle_type_XXX
---------------------------------------------

Each particle type is specified in another subgroup `particle_type_XXX`, see :ref:`particle_model_config`.


Group /input/model/unit_XXX/discretization - UNIT_TYPE - COLUMN_MODEL_1D
-------------------------------------------------------------------------

``USE_ANALYTIC_JACOBIAN``

   Determines whether analytically computed Jacobian matrix (faster) is used (value is 1) instead of Jacobians generated by algorithmic differentiation (slower, value is 0)
   
   =============  ===========================  =============
   **Type:** int  **Range:** :math:`\{0, 1\}`  **Length:** 1
   =============  ===========================  =============

Spatial discretization - Numerical Methods
------------------------------------------

CADET offers two spatial discretization methods: Finite Volumes (FV) and Discontinuous Galerkin (DG). Each method has it's own set of input fields.
While both methods approximate the same solution to the same underlying model, they may differ in terms of computational performance.
With our currently implemented variants of FV and DG, FV perform better for solutions with steep gradients or discontinuities, while DG can be much faster for rather smooth solutions.
For the same number of discrete points, DG will generally be slower but often more accurate.

For further information on the choice of discretization methods and their parameters, see :ref:`spatial_discretization_methods`.

``SPATIAL_METHOD``

   Spatial discretization method. Optional, defaults to :math:`\texttt{FV}`

   ================  ===============================================  =============
   **Type:** string  **Range:** :math:`\{\texttt{FV}, \texttt{DG}\}`  **Length:** 1
   ================  ===============================================  =============

Discontinuous Galerkin
----------------------

``POLYDEG``

   DG polynomial degree. Optional, defaults to 4 and :math:`N_d \in \{3, 4, 5\}` is recommended. The total number of axial discrete points is given by (``POLYDEG`` + 1 ) * ``NELEM``
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``NELEM``

   Number of axial column discretization DG cells\elements. The total number of axial discrete points is given by (``POLYDEG`` + 1 ) * ``NELEM``
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``NCOL``

   Number of axial discrete points. Optional and ignored if ``NELEM`` is defined. Otherwise, used to calculate ``NELEM`` = :math:`\lfloor` ``NCOL`` / (``POLYDEG`` + 1 ) :math:`\rfloor`
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``EXACT_INTEGRATION``

   Specifies the DG integration variant. Optional, defaults to 0
   
   =============  ===========================  =============
   **Type:** int  **Range:** :math:`\{0, 1\}`  **Length:** 1
   =============  ===========================  =============

``LINEAR_SOLVER``

   Specifies the linear solver variant used to factorize the semidiscretized system. Optional, defaults to ``SparseLU``. For more information on these solvers, we refer to the `Eigen documentation <https://eigen.tuxfamily.org/>`_
   
   =============  ===================================================================================  =============
   **Type:** int  **Range:** :math:`\{\texttt{SparseLU}, \texttt{SparseQR}, ..., \texttt{BiCGSTAB}\}`  **Length:** 1
   =============  ===================================================================================  =============
   

   When using the DG method, we generally recommend specifying ``USE_MODIFIED_NEWTON = 1`` in :ref:`FFSolverTime`, i.e. to use the modified Newton method to solve the linear system within the time integrator.
   For further information on discretization parameters, see also :ref:`non_consistency_solver_parameters`.

Finite Volumes
--------------

``NCELLS``

   Number of axial column discretization points, i.e. FV cells
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``RECONSTRUCTION``

   Type of reconstruction method for fluxes
   
   ================  ================================  =============
   **Type:** string  **Range:** :math:`\texttt{WENO}`  **Length:** 1
   ================  ================================  =============

The following FV discretization parameters are only required if particles are present:

``GS_TYPE``

   Type of Gram-Schmidt orthogonalization, see IDAS guide Section 4.5.7.3, p. 41f. A value of :math:`0` enables classical Gram-Schmidt, a value of 1 uses modified Gram-Schmidt.
   
   =============  ===========================  =============
   **Type:** int  **Range:** :math:`\{0, 1\}`  **Length:** 1
   =============  ===========================  =============

``MAX_KRYLOV``

   Defines the size of the Krylov subspace in the iterative linear GMRES solver (0: :math:`\texttt{MAX_KRYLOV} = \texttt{NCOL} \cdot \texttt{NCOMP} \cdot \texttt{NPARTYPE}`)
   
   =============  ============================================================================================  =============
   **Type:** int  **Range:** :math:`\{0, \dots, \texttt{NCOL} \cdot \texttt{NCOMP} \cdot \texttt{NPARTYPE} \}`  **Length:** 1
   =============  ============================================================================================  =============

``MAX_RESTARTS``

   Maximum number of restarts in the GMRES algorithm. If lack of memory is not an issue, better use a larger Krylov space than restarts.
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 0`  **Length:** 1
   =============  =========================  =============

``SCHUR_SAFETY``

   Schur safety factor; Influences the tradeoff between linear iterations and nonlinear error control; see IDAS guide Section~2.1 and 5.
   
   ================  =========================  =============
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** 1
   ================  =========================  =============

``FIX_ZERO_SURFACE_DIFFUSION``

   Determines whether the surface diffusion parameters :math:`\texttt{SURFACE_DIFFUSION}` are fixed if the parameters are zero. If the parameters are fixed to zero (:math:`\texttt{FIX_ZERO_SURFACE_DIFFUSION} = 1`, :math:`\texttt{SURFACE_DIFFUSION} = 0`), the parameters must not become non-zero during this or subsequent simulation runs. The internal data structures are optimized for a more efficient simulation.  This field is optional and defaults to :math:`0` (optimization disabled in favor of flexibility).
   
   =============  ===========================  =============
   **Type:** int  **Range:** :math:`\{0, 1\}`  **Length:** 1
   =============  ===========================  =============

For further information on discretization parameters, see also :ref:`flux_reconstruction_methods` (FV specific)), and :ref:`non_consistency_solver_parameters`.
