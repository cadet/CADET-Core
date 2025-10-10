.. _axial_flow_column_2D_config:

Axial Flow Column 2D
====================

Group /input/model/unit_XXX - UNIT_TYPE - COLUMN_MODEL_2D
---------------------------------------------------------

``UNIT_TYPE``

   Specifies the type of unit operation model
   
   ================  ===========================================  =============
   **Type:** string  **Range:** :math:`\texttt{COLUMN_MODEL_2D}`  **Length:** 1
   ================  ===========================================  =============

``NCOMP``

   Number of chemical components in the chromatographic medium
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``CROSS_SECTION_AREA``

   Cross section area of the column. This parameter is optional and will be ignored if `COL_RADIUS` is provided.
   **Unit:** :math:`\mathrm{m}^{2}`
   
   ================  =====================  =============
   **Type:** double  **Range:** :math:`>0`  **Length:** 1
   ================  =====================  =============

``COL_LENGTH``

   Column length / height

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

``COL_POROSITY``

   Column porosity, either constant (length is 1) or for each radial zone (length is :math:`\texttt{NRAD}`). 
   In case of a spatially inhomogeneous setting, the :math:`\texttt{SENS_PARTYPE}` field is used for indexing the radial zone when specifying parameter sensitivities.

   ================  ========================  =====================================
   **Type:** double  **Range:** :math:`(0,1]`  **Length:** :math:`1 / \texttt{NRAD}`
   ================  ========================  =====================================

``PAR_TYPE_VOLFRAC``

   Volume fractions of the particle types. The volume fractions can be set homogeneous or individually along both axes. For each cell, the volume fractions have to sum to 1. 
   In case of a spatially inhomogeneous setting, the :math:`\texttt{SENS_SECTION}` field is used for indexing the axial cell and the :math:`\texttt{SENS_REACTION}` field is used for indexing the radial cell when specifying parameter sensitivities.  This field is optional in case of only one particle type.

   ================  ========================  ===========================================================
   **Type:** double  **Range:** :math:`[0,1]`  **Length:** see :math:`\texttt{PAR_TYPE_VOLFRAC_MULTIPLEX}`
   ================  ========================  ===========================================================

``PAR_TYPE_VOLFRAC_MULTIPLEX``

   Multiplexing mode of :math:`\texttt{PAR_TYPE_VOLFRAC}`. Determines whether :math:`\texttt{PAR_TYPE_VOLFRAC}` is treated as radial- and/or section-independent.
   This field is optional. When left out, multiplexing behavior is inferred from the length of :math:`\texttt{PAR_TYPE_VOLFRAC}`.  Valid modes are:

  0. Radial-independent, axial-independent; length of :math:`\texttt{PAR_TYPE_VOLFRAC}` is :math:`\texttt{NPARTYPE}`
  1. Radial-dependent, axial-independent; length of :math:`\texttt{PAR_TYPE_VOLFRAC}` is :math:`\texttt{NRAD} \cdot \texttt{NPARTYPE}`; ordering is radial-major
  2. Axial-dependent; length of :math:`\texttt{PAR_TYPE_VOLFRAC}` is :math:`\texttt{NCOL} \cdot \texttt{NPARTYPE}`; ordering is axial-major
  3. Radial-dependent, axial-dependent; length of :math:`\texttt{PAR_TYPE_VOLFRAC}` is :math:`\texttt{NCOL} \cdot \texttt{NRAD} \cdot \texttt{NPARTYPE}`; ordering is axial-radial-major

``VELOCITY``

   Indicates flow direction in each radial zone (forward if value is positive, backward if value is negative), see Section :ref:`MUOPGRMflow2D`).  In case of a spatially inhomogeneous setting, the :math:`\texttt{SENS_PARTYPE}` field is used for indexing the radial cell when specifying parameter sensitivities.

   ================  =============================  ===================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** see :math:`\texttt{VELOCITY_MULTIPLEX}`
   ================  =============================  ===================================================

``VELOCITY_MULTIPLEX``

   Multiplexing mode of :math:`\texttt{VELOCITY}`. Determines whether :math:`\texttt{VELOCITY}` is treated as radial- and/or section-independent.  This field is optional. When left out, multiplexing behavior is inferred from the length of :math:`\texttt{VELOCITY}`.  Valid modes are:

  0. Radial-independent, section-independent; length of :math:`\texttt{VELOCITY}` is 1
  1. Radial-dependent, section-independent; length of :math:`\texttt{VELOCITY}` is :math:`\texttt{NRAD}`
  2. Section-dependent; length of :math:`\texttt{VELOCITY}` is :math:`\texttt{NSEC}`
  3. Radial-dependent, section-dependent; length of :math:`\texttt{VELOCITY}` is :math:`\texttt{NRAD} \cdot \texttt{NSEC}`; ordering is section-major

``COL_DISPERSION_AXIAL``

   Axial dispersion coefficient

   **Unit:** :math:`\mathrm{m}_{\mathrm{IV}}^{2}\,\mathrm{s}^{-1}`
   
   ================  =========================  ===============================================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** see :math:`\texttt{COL_DISPERSION_AXIAL_MULTIPLEX}`
   ================  =========================  ===============================================================

``COL_DISPERSION_AXIAL_MULTIPLEX``

   Multiplexing mode of :math:`\texttt{COL_DISPERSION_AXIAL}`. Determines whether :math:`\texttt{COL_DISPERSION_AXIAL}` is treated as component-, radial-, and/or section-independent.  This field is optional. When left out, multiplexing behavior is inferred from the length of :math:`\texttt{COL_DISPERSION_AXIAL}`.  Valid modes are:

  0. Component-independent, radial-independent, section-independent; length of :math:`\texttt{COL_DISPERSION_AXIAL}` is 1
  1. Component-independent, radial-dependent, section-independent; length of :math:`\texttt{COL_DISPERSION_AXIAL}` is :math:`\texttt{NRAD}`
  2. Component-dependent, radial-independent, section-independent; length of :math:`\texttt{COL_DISPERSION_AXIAL}` is :math:`\texttt{NCOMP}`
  3. Component-dependent, radial-dependent, section-independent; length of :math:`\texttt{COL_DISPERSION_AXIAL}` is :math:`\texttt{NCOMP} \cdot \texttt{NRAD}`; ordering is radial-major
  4. Component-independent, radial-independent, section-dependent; length of :math:`\texttt{COL_DISPERSION_AXIAL}` is :math:`\texttt{NSEC}`
  5. Component-independent, radial-dependent, section-dependent; length of :math:`\texttt{COL_DISPERSION_AXIAL}` is :math:`\texttt{NRAD} \cdot \texttt{NSEC}`; ordering is section-major
  6. Component-dependent, radial-independent, section-independent; length of :math:`\texttt{COL_DISPERSION_AXIAL}` is :math:`\texttt{NCOMP} \cdot \texttt{NSEC}`; ordering is section-major
  7. Component-dependent, radial-dependent, section-dependent; length of :math:`\texttt{COL_DISPERSION_AXIAL}` is :math:`\texttt{NCOMP} \cdot \texttt{NRAD} \cdot \texttt{NSEC}`; ordering is section-radial-major
   
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
  1. Component-independent, radial-dependent, section-independent; length of :math:`\texttt{COL_DISPERSION_RADIAL}` is :math:`\texttt{NRAD}`
  2. Component-dependent, radial-independent, section-independent; length of :math:`\texttt{COL_DISPERSION_RADIAL}` is :math:`\texttt{NCOMP}`
  3. Component-dependent, radial-dependent, section-independent; length of :math:`\texttt{COL_DISPERSION_RADIAL}` is :math:`\texttt{NCOMP} \cdot \texttt{NRAD}`; ordering is radial-major
  4. Component-independent, radial-independent, section-dependent; length of :math:`\texttt{COL_DISPERSION_RADIAL}` is :math:`\texttt{NSEC}`
  5. Component-independent, radial-dependent, section-dependent; length of :math:`\texttt{COL_DISPERSION_RADIAL}` is :math:`\texttt{NRAD} \cdot \texttt{NSEC}`; ordering is section-major
  6. Component-dependent, radial-independent, section-independent; length of :math:`\texttt{COL_DISPERSION_RADIAL}` is :math:`\texttt{NCOMP} \cdot \texttt{NSEC}`; ordering is section-major
  7. Component-dependent, radial-dependent, section-dependent; length of :math:`\texttt{COL_DISPERSION_RADIAL}` is :math:`\texttt{NCOMP} \cdot \texttt{NRAD} \cdot \texttt{NSEC}`; ordering is section-radial-major

   =============  ===================================  =============
   **Type:** int  **Range:** :math:`\{0, \dots, 7 \}`  **Length:** 1
   =============  ===================================  =============

``REACTION_MODEL_BULK``

   Specifies the type of reaction model of the bulk volume. The model is configured in the subgroup :math:`\texttt{reaction_bulk}`.
   
   ================  ========================================  =============
   **Type:** string  **Range:** See Section :ref:`FFReaction`  **Length:** 1
   ================  ========================================  =============

``INIT_C``

   Initial concentrations for each component in all radial zones the bulk mobile phase (length :math:`\texttt{NCOMP}`), or for each component in each radial zone (length :math:`\texttt{NCOMP} \cdot \texttt{NRAD}`, ordering radial-major)

   **Unit:** :math:`\mathrm{mol}\,\mathrm{m}_{\mathrm{IV}}^{-3}`

   ================  =========================  =========================================================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NCOMP} / \texttt{NCOMP} \cdot \texttt{NRAD}`
   ================  =========================  =========================================================================

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


Group /input/model/unit_XXX/discretization - UNIT_TYPE - COLUMN_MODEL_2D
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

``RADIAL_DISC_TYPE``

   Specifies the radial discretization scheme. Valid values are :math:`\texttt{EQUIDISTANT}`, :math:`\texttt{EQUIVOLUME}`, and :math:`\texttt{USER_DEFINED}`.

   ================  =============
   **Type:** string  **Length:** 1
   ================  =============

``RADIAL_COMPARTMENTS``

   Coordinates for the radial compartment boundaries (ignored if :math:`\texttt{RADIAL_DISC_TYPE} \neq \texttt{USER_DEFINED}`). The coordinates are absolute and have to include the endpoints 0 and :math:`\texttt{COLUMN_RADIUS}`. The values are expected in ascending order (i.e., 0 is the first and :math:`\texttt{COLUMN_RADIUS}` the last value in the array).

   **Unit:** :math:`\mathrm{m}`

   ================  =============================================  ====================================
   **Type:** double  **Range:** :math:`[0,\texttt{COLUMN_RADIUS}]`  **Length:** :math:`\texttt{NRAD} + 1`
   ================  =============================================  ====================================

Discontinuous Galerkin
----------------------

``AX_POLYDEG``

   DG polynomial degree for axial discretization. Optional, defaults to 4 and :math:`N^z_d \in \{3, 4, 5\}` is recommended.
   The total number of axial discrete points is given by (``AX_POLYDEG`` + 1 ) * ``AX_NELEM``
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``AX_NELEM``

   Number of axial column discretization DG cells\elements. The total number of axial discrete points is given by (``POLYDEG`` + 1 ) * ``NELEM``
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``RAD_POLYDEG``

   DG polynomial degree for radial discretization. Optional, defaults to 4 and :math:`N^r_d \in \{3, 4, 5\}` is recommended, and should generally be the same as the axial degree.
   The total number of radial discrete points is given by (``RAD_POLYDEG`` + 1 ) * ``RAD_NELEM``
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``RAD_NELEM``

   Number of radial column discretization DG cells\elements. The total number of axial discrete points is given by (``POLYDEG`` + 1 ) * ``NELEM``
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``LINEAR_SOLVER``

   Specifies the linear solver variant used to factorize the semidiscretized system. Optional, defaults to ``SparseLU``. For more information on these solvers, we refer to the `Eigen documentation <https://eigen.tuxfamily.org/>`_
   
   =============  ===================================================================================  =============
   **Type:** int  **Range:** :math:`\{\texttt{SparseLU}, \texttt{SparseQR}, ..., \texttt{BiCGSTAB}\}`  **Length:** 1
   =============  ===================================================================================  =============
   

   For further information on discretization parameters, see also :ref:`non_consistency_solver_parameters`.

Finite Volumes
--------------

``NCOL``

   Number of axial column discretization cells

   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``NRAD``

   Number of radial column discretization cells

   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============

``NPAR``

   Number of particle (radial) discretization cells for each particle type

   =============  =========================  =========================================
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** :math:`1 / \texttt{NPARTYPE}`
   =============  =========================  =========================================

``LINEAR_SOLVER_BULK``

   Linear solver used for the sparse column bulk block. This field is optional, the best available method is selected (i.e., sparse direct solver if possible).  Valid values are:

  - :math:`\texttt{DENSE}` Converts the sparse matrix into a banded matrix and uses regular LAPACK. Slow and memory intensive, but always available.
  - :math:`\texttt{UMFPACK}` Uses the UMFPACK sparse direct solver (LU decomposition) from SuiteSparse. Fast, but has to be enabled when compiling and requires UMFPACK library.
  - :math:`\texttt{SUPERLU}` Uses the SuperLU sparse direct solver (LU decomposition). Fast, but has to be enabled when compiling and requires SuperLU library.

   ================  =======================================================================  =============
   **Type:** string  **Range:** :math:`\{\texttt{DENSE},\texttt{UMFPACK},\texttt{SUPERLU}\}`  **Length:** 1
   ================  =======================================================================  =============

``RECONSTRUCTION``

   Type of reconstruction method for fluxes

   ================  ================================  =============
   **Type:** string  **Range:** :math:`\texttt{WENO}`  **Length:** 1
   ================  ================================  =============

``GS_TYPE``

   Type of Gram-Schmidt orthogonalization, see IDAS guide Section~4.5.7.3, p.~41f. A value of 0 enables classical Gram-Schmidt, a value of 1 uses modified Gram-Schmidt.

   =============  ===========================  =============
   **Type:** int  **Range:** :math:`\{0, 1\}`  **Length:** 1
   =============  ===========================  =============

``MAX_KRYLOV``

   Defines the size of the Krylov subspace in the iterative linear GMRES solver (0: :math:`\texttt{MAX_KRYLOV} = \texttt{NCOL} \cdot \texttt{NRAD} \cdot \texttt{NCOMP} \cdot \texttt{NPARTYPE}`)

   =============  ================================================================================================================  =============
   **Type:** int  **Range:** :math:`\{0, \dots, \texttt{NCOL} \cdot \texttt{NRAD} \cdot \texttt{NCOMP} \cdot \texttt{NPARTYPE} \}`  **Length:** 1
   =============  ================================================================================================================  =============

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

For further discretization parameters, see also :ref:`flux_reconstruction_methods`, and :ref:`non_consistency_solver_parameters`.
