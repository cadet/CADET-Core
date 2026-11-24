(radial-flow-column-1d-config)=
# Radial Flow Column 1D

## Group /input/model/unit_XXX - UNIT_TYPE - COLUMN_MODEL_1D

`UNIT_TYPE`

: Specifies the type of unit operation model

  ```{list-table}

  * - **Type:** string
    - **Range:** $\texttt{COLUMN_MODEL_1D}$
    - **Length:** 1
  ```

`GEOMETRY`

: Column geometry

  ```{list-table}

  * - **Type:** string
    - **Range:** $\texttt{RADIAL_FLOW_CYLINDER_SHELL}$
    - **Length:** 1
  ```

`NCOMP`

: Number of chemical components in the chromatographic medium

  ```{list-table}

  * - **Type:** int
    - **Range:** $\geq 1$
    - **Length:** 1
  ```

`CROSS_SECTION_AREA_INNER`

: Cross section area of the inner column radius (optional, only triggers a geometry check against `CROSS_SECTION_AREA_OUTER`, `CYLINDER_HEIGHT`, and `BED_LENGTH`).
  The inner radius has to be positive, i.e. `BED_LENGTH` strictly smaller than the outer radius:
  the interstitial velocity of a radial flow column scales with $1/r$ and is singular on the axis.

  **Unit:** $\mathrm{m}^2$

  ```{list-table}

  * - **Type:** double
    - **Range:** $> 0$
    - **Length:** 1
  ```

`CROSS_SECTION_AREA_OUTER`

: Cross section area of the outer column radius

  **Unit:** $\mathrm{m}^2$

  ```{list-table}

  * - **Type:** double
    - **Range:** $> 0$
    - **Length:** 1
  ```

`BED_LENGTH`

: Column length/height, i.e., outer minus inner column radius.

  **Unit:** $\mathrm{m}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $> 0$
    - **Length:** 1
  ```

`CYLINDER_HEIGHT`

: Height of the cylinder. Note: not the bed length

  **Unit:** $\mathrm{m}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $> 0$
    - **Length:** 1
  ```

`FORWARD_FLOW`

: Flow direction of the column. If set to 1, the flow is from the outer radius to the inner radius.

  ```{list-table}

  * - **Type:** bool
    - **Range:** $\{0,1\}$
    - **Length:** 1 or NSEC
  ```

`CIRCLE_FRACTION`

: Optional, only required for wedge geometry. Specifies the cylinder shell fraction of the wedge

  ```{list-table}

  * - **Type:** double
    - **Range:** $(0,1]$
    - **Length:** 1
  ```

`COL_POROSITY`

: Column porosity

  ```{list-table}

  * - **Type:** double
    - **Range:** $(0,1]$
    - **Length:** 1
  ```

`NPARTYPE`

: Number of particle types.

  ```{list-table}

  * - **Type:** int
    - **Range:** $\geq 0$
    - **Length:** 1
  ```

`PAR_TYPE_VOLFRAC`

: Volume fractions of the particle types. The volume fractions can be set for all axial cells together or for each individual axial cell. For each cell, the volume fractions have to sum to $1$. In case of a spatially inhomogeneous setting, the data is expected in cell-major ordering and the $\texttt{SENS_SECTION}$ field is used for indexing the axial cell when specifying parameter sensitivities. This field is optional in case of only one particle type.

  ```{list-table}

  * - **Type:** double
    - **Range:** $[0,1]$
    - **Length:** $\texttt{NPARTYPE} / \texttt{NCOL} \cdot \texttt{NPARTYPE}$
  ```

`COL_DISPERSION`

: Axial dispersion coefficient

  **Unit:** $\mathrm{m}_{\mathrm{IV}}^{2}\,\mathrm{s}^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\geq 0$
    - **Length:** see $\texttt{COL_DISPERSION_MULTIPLEX}$
  ```

`COL_DISPERSION_MULTIPLEX`

: Multiplexing mode of $\texttt{COL_DISPERSION}$. Determines whether $\texttt{COL_DISPERSION}$ is treated as component- and/or section-independent. This field is optional. When left out, multiplexing behavior is inferred from the length of $\texttt{COL_DISPERSION}$. Valid modes are:

  0. Component-independent, section-independent; length of $\texttt{COL_DISPERSION}$ is $1$
  1. Component-dependent, section-independent; length of $\texttt{COL_DISPERSION}$ is $\texttt{NCOMP}$
  2. Component-independent, section-dependent; length of $\texttt{COL_DISPERSION}$ is $\texttt{NSEC}$
  3. Component-dependent, section-dependent; length of $\texttt{COL_DISPERSION}$ is $\texttt{NCOMP} \cdot \texttt{NSEC}$; ordering is section-major

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0, \dots, 3 \}$
    - **Length:** 1
  ```

`NREAC_LIQUID`

: Number of liquid phase reaction models (optional, only if liquid reactions are present).

  ```{list-table}

  * - **Type:** int
    - **Range:** $\geq 0$
    - **Length:** 1
  ```

`INIT_C`

: Initial concentrations for each component in the bulk mobile phase

  **Unit:** $\mathrm{mol}\,\mathrm{m}_{\mathrm{IV}}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\geq 0$
    - **Length:** $\texttt{NCOMP}$
  ```

`INIT_STATE`

: Full state vector for initialization (optional, $\texttt{INIT_C}$, $\texttt{INIT_CP}$, and $\texttt{INIT_CS}$ will be ignored; if length is $2\texttt{NDOF}$, then the second half is used for time derivatives).
  The ordering of the state vector is defined in {ref}`unitoperationstateordering`.

  **Unit:** $various$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** $\texttt{NDOF} / 2\texttt{NDOF}$
  ```

## Group /input/model/unit_XXX/liquid_reaction_YYY

Each reaction is specified in another subgroup `liquid_reaction_YYY`, see {ref}`ffreaction`.

## Group /input/model/unit_XXX/particle_type_XXX

Each particle type is specified in another subgroup `particle_type_XXX`, see {ref}`particle-model-config`.

## Group /input/model/unit_XXX/discretization - UNIT_TYPE - COLUMN_MODEL_1D

`USE_ANALYTIC_JACOBIAN`

: Determines whether analytically computed Jacobian matrix (faster) is used (value is 1) instead of Jacobians generated by algorithmic differentiation (slower, value is 0)

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0, 1\}$
    - **Length:** 1
  ```

## Spatial discretization - Numerical Methods

CADET offers two spatial discretization methods: Finite Volumes (FV) and Discontinuous Galerkin (DG). Each method has its own set of input fields.
While both methods approximate the same solution to the same underlying model, they may differ in terms of computational performance.

For further information on the choice of discretization methods and their parameters, see {ref}`spatial-discretization-methods`.

`SPATIAL_METHOD`

: Spatial discretization method

  ```{list-table}

  * - **Type:** string
    - **Range:** $\{\texttt{FV}, \texttt{DG}\}$
    - **Length:** 1
  ```

### Discontinuous Galerkin

`POLYDEG`

: DG polynomial degree. Optional, defaults to 4 and $N_d \in \{3, 4, 5\}$ is recommended. The total number of radial discrete points is given by (`POLYDEG` + 1 ) * `NELEM`

  ```{list-table}

  * - **Type:** int
    - **Range:** $\geq 1$
    - **Length:** 1
  ```

`NELEM`

: Number of radial column discretization DG cells/elements. The total number of radial discrete points is given by (`POLYDEG` + 1 ) * `NELEM`

  ```{list-table}

  * - **Type:** int
    - **Range:** $\geq 1$
    - **Length:** 1
  ```

`NCOL`

: Number of radial discrete points. Optional and ignored if `NELEM` is defined. Otherwise, used to calculate `NELEM` = $\lfloor$ `NCOL` / (`POLYDEG` + 1 ) $\rfloor$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\geq 1$
    - **Length:** 1
  ```

`POLYNOMIAL_INTERPOLATION_NODES`

: Specifies the interpolation nodes for the DG ansatz. Optional, defaults to `LEGENDRE_GAUSS_LOBATTO`.

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{\texttt{LEGENDRE_GAUSS_LOBATTO}, \texttt{CHEBYSHEV_GAUSS_LOBATTO}\}$
    - **Length:** 1
  ```

`LINEAR_SOLVER`

: Specifies the linear solver variant used to factorize the semidiscretized system. Optional, defaults to `SparseLU`. For more information on these solvers, we refer to the [Eigen documentation](https://eigen.tuxfamily.org/)

  ```{list-table}

  * - **Type:** string
    - **Range:** $\{\texttt{SparseLU}, \texttt{SparseQR}, ..., \texttt{BiCGSTAB}\}$
    - **Length:** 1
  ```

When using the DG method, we generally recommend specifying `USE_MODIFIED_NEWTON = 1` in {ref}`ffsolvertime`, i.e. to use the modified Newton method to solve the linear system within the time integrator.

### Finite Volumes

`NCOL`

: Number of radial column discretization points, i.e. FV cells

  ```{list-table}

  * - **Type:** int
    - **Range:** $\geq 1$
    - **Length:** 1
  ```

`RECONSTRUCTION`

: Type of reconstruction method for fluxes

  ```{list-table}

  * - **Type:** string
    - **Range:** $\texttt{WENO, KOREN}$
    - **Length:** 1
  ```

  For further information on discretization parameters for reconstruction methods, see also {ref}`flux-reconstruction-methods` (FV specific).

The following FV discretization parameters configure the GMRES solver of the unit operation's own
arrow-head Schur complement. They are required if particles are present **and**
`FV_ARROW_HEAD_OPTIMIZATION` $= 1$, which is the default for an FV bulk with FV particles.
With `FV_ARROW_HEAD_OPTIMIZATION` $= 0$ the unit is solved by a global sparse solver that
does not read them. They are independent of the identically named fields in
{ref}`ffmodelsolver`, which size the GMRES of the network.

`GS_TYPE`

: Type of Gram-Schmidt orthogonalization, see IDAS guide Section 4.5.7.3, p. 41f. A value of $0$ enables classical Gram-Schmidt, a value of 1 uses modified Gram-Schmidt.

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0, 1\}$
    - **Length:** 1
  ```

`MAX_KRYLOV`

: Defines the size of the Krylov subspace in the iterative linear GMRES solver (0: $\texttt{MAX_KRYLOV} = \texttt{NCOL} \cdot \texttt{NCOMP} \cdot \texttt{NPARTYPE}$)

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0, \dots, \texttt{NCOL} \cdot \texttt{NCOMP} \cdot \texttt{NPARTYPE} \}$
    - **Length:** 1
  ```

`MAX_RESTARTS`

: Maximum number of restarts in the GMRES algorithm. If lack of memory is not an issue, better use a larger Krylov space than restarts.

  ```{list-table}

  * - **Type:** int
    - **Range:** $\geq 0$
    - **Length:** 1
  ```

`SCHUR_SAFETY`

: Schur safety factor; Influences the tradeoff between linear iterations and nonlinear error control; see IDAS guide Section 2.1 and 5.

  ```{list-table}

  * - **Type:** double
    - **Range:** $\geq 0$
    - **Length:** 1
  ```

`FV_ARROW_HEAD_OPTIMIZATION`

: Optional, specifies whether an optimized Jacobian block structure is exploited by a customized linear solver module for improved computational efficiency.
  Prohibits mix of particle types and mix of FV, DG discretizations for bulk and particles.

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
    - **Length:** 1
  ```

When using the FV method, we generally recommend specifying `USE_MODIFIED_NEWTON = 0` in {ref}`ffsolvertime`, i.e. to use the full Newton method to solve the linear system within the time integrator.
For further information on discretization parameters, see also {ref}`flux-reconstruction-methods` (FV specific), and {ref}`non-consistency-solver-parameters`.
