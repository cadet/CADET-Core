(axial-flow-column-2d-config)=

# Axial Flow Column 2D

## Group /input/model/unit_XXX - UNIT_TYPE - COLUMN_MODEL_2D

`UNIT_TYPE`

> Specifies the type of unit operation model
>
> ```{list-table}
>
> * - **Type:** string
>   - **Range:** $\texttt{COLUMN_MODEL_2D}$
>   - **Length:** 1
> ```

`NCOMP`

> Number of chemical components in the chromatographic medium
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq 1$
>   - **Length:** 1
> ```

`CROSS_SECTION_AREA`

> Cross section area of the column. This parameter is only read by the FV discretization and is ignored if `COL_RADIUS` is provided.
>
> **Unit:** $\mathrm{m}^{2}$
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $>0$
>   - **Length:** 1
> ```

`BED_LENGTH`

> Column length / height
>
> **Unit:** $\mathrm{m}$
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $> 0$
>   - **Length:** 1
> ```

`COL_RADIUS`

> Column radius. Required for the DG discretization; for the FV discretization it is optional if `CROSS_SECTION_AREA` is provided.
>
> **Unit:** $\mathrm{m}$
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $> 0$
>   - **Length:** 1
> ```

`COL_POROSITY`

> Column porosity, either constant (length is 1) or for each radial zone (length is $\texttt{RAD_NELEM}$).
> In case of a spatially inhomogeneous setting, the $\texttt{SENS_PARTYPE}$ field is used for indexing the radial zone when specifying parameter sensitivities.
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $(0,1]$
>   - **Length:** $1 / \texttt{RAD_NELEM}$
> ```

`NPARTYPE`

> Number of particle types. A value of $0$ means that no particles are present.
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq 0$
>   - **Length:** 1
> ```

`PAR_TYPE_VOLFRAC`

> Volume fractions of the particle types. The volume fractions can be set homogeneous or individually along both axes. For each cell, the volume fractions have to sum to 1.
> In case of a spatially inhomogeneous setting, the $\texttt{SENS_SECTION}$ field is used for indexing the axial cell and the $\texttt{SENS_REACTION}$ field is used for indexing the radial cell when specifying parameter sensitivities. This field is optional in case of only one particle type.
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $[0,1]$
>   - **Length:** see $\texttt{PAR_TYPE_VOLFRAC_MULTIPLEX}$
> ```

`PAR_TYPE_VOLFRAC_MULTIPLEX`

> Multiplexing mode of $\texttt{PAR_TYPE_VOLFRAC}$. Determines whether $\texttt{PAR_TYPE_VOLFRAC}$ is treated as radial- and/or section-independent.
> This field is optional. When left out, multiplexing behavior is inferred from the length of $\texttt{PAR_TYPE_VOLFRAC}$. Valid modes are:
>
> 0. Radial-independent, axial-independent; length of $\texttt{PAR_TYPE_VOLFRAC}$ is $\texttt{NPARTYPE}$
> 1. Radial-dependent, axial-independent; length of $\texttt{PAR_TYPE_VOLFRAC}$ is $\texttt{RAD_NELEM} \cdot \texttt{NPARTYPE}$; ordering is radial-major
> 2. Axial-dependent; length of $\texttt{PAR_TYPE_VOLFRAC}$ is $\texttt{NCOL} \cdot \texttt{NPARTYPE}$; ordering is axial-major
> 3. Radial-dependent, axial-dependent; length of $\texttt{PAR_TYPE_VOLFRAC}$ is $\texttt{NCOL} \cdot \texttt{RAD_NELEM} \cdot \texttt{NPARTYPE}$; ordering is axial-radial-major
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\{0, \dots, 3 \}$
>   - **Length:** 1
> ```

`VELOCITY`

> Interstitial velocity and flow direction in each radial zone (forward if value is positive, backward if value is negative), see Section {ref}`muopgrmflow2d`. Optional; if left out, the velocity is inferred from the volumetric flow rates of the network connections. In case of a spatially inhomogeneous setting, the $\texttt{SENS_PARTYPE}$ field is used for indexing the radial cell when specifying parameter sensitivities.
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $\mathbb{R}$
>   - **Length:** see $\texttt{VELOCITY_MULTIPLEX}$
> ```

`VELOCITY_MULTIPLEX`

> Multiplexing mode of $\texttt{VELOCITY}$. Determines whether $\texttt{VELOCITY}$ is treated as radial- and/or section-independent. This field is optional. When left out, multiplexing behavior is inferred from the length of $\texttt{VELOCITY}$. Valid modes are:
>
> 0. Radial-independent, section-independent; length of $\texttt{VELOCITY}$ is 1
> 1. Radial-dependent, section-independent; length of $\texttt{VELOCITY}$ is $\texttt{RAD_NELEM}$
> 2. Section-dependent; length of $\texttt{VELOCITY}$ is $\texttt{NSEC}$
> 3. Radial-dependent, section-dependent; length of $\texttt{VELOCITY}$ is $\texttt{RAD_NELEM} \cdot \texttt{NSEC}$; ordering is section-major
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\{0, \dots, 3 \}$
>   - **Length:** 1
> ```

`COL_DISPERSION_AXIAL`

> Axial dispersion coefficient
>
> **Unit:** $\mathrm{m}_{\mathrm{IV}}^{2}\,\mathrm{s}^{-1}$
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $\geq 0$
>   - **Length:** see $\texttt{COL_DISPERSION_AXIAL_MULTIPLEX}$
> ```

`COL_DISPERSION_AXIAL_MULTIPLEX`

> Multiplexing mode of $\texttt{COL_DISPERSION_AXIAL}$. Determines whether $\texttt{COL_DISPERSION_AXIAL}$ is treated as component-, radial-, and/or section-independent. This field is optional. When left out, multiplexing behavior is inferred from the length of $\texttt{COL_DISPERSION_AXIAL}$. Valid modes are:
>
> 0. Component-independent, radial-independent, section-independent; length of $\texttt{COL_DISPERSION_AXIAL}$ is 1
> 1. Component-independent, radial-dependent, section-independent; length of $\texttt{COL_DISPERSION_AXIAL}$ is $\texttt{RAD_NELEM}$
> 2. Component-dependent, radial-independent, section-independent; length of $\texttt{COL_DISPERSION_AXIAL}$ is $\texttt{NCOMP}$
> 3. Component-dependent, radial-dependent, section-independent; length of $\texttt{COL_DISPERSION_AXIAL}$ is $\texttt{NCOMP} \cdot \texttt{RAD_NELEM}$; ordering is radial-major
> 4. Component-independent, radial-independent, section-dependent; length of $\texttt{COL_DISPERSION_AXIAL}$ is $\texttt{NSEC}$
> 5. Component-independent, radial-dependent, section-dependent; length of $\texttt{COL_DISPERSION_AXIAL}$ is $\texttt{RAD_NELEM} \cdot \texttt{NSEC}$; ordering is section-major
> 6. Component-dependent, radial-independent, section-dependent; length of $\texttt{COL_DISPERSION_AXIAL}$ is $\texttt{NCOMP} \cdot \texttt{NSEC}$; ordering is section-major
> 7. Component-dependent, radial-dependent, section-dependent; length of $\texttt{COL_DISPERSION_AXIAL}$ is $\texttt{NCOMP} \cdot \texttt{RAD_NELEM} \cdot \texttt{NSEC}$; ordering is section-radial-major
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\{0, \dots, 7 \}$
>   - **Length:** 1
> ```

`COL_DISPERSION_RADIAL`

> Radial dispersion coefficient. In case of a spatially inhomogeneous setting, the $\texttt{SENS_PARTYPE}$ field is used for indexing the radial zone when specifying parameter sensitivities.
>
> **Unit:** $\mathrm{m}_{\mathrm{IV}}^{2}\,\mathrm{s}^{-1}$
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $\geq 0$
>   - **Length:** see $\texttt{COL_DISPERSION_RADIAL_MULTIPLEX}$
> ```

`COL_DISPERSION_RADIAL_MULTIPLEX`

> Multiplexing mode of $\texttt{COL_DISPERSION_RADIAL}$. Determines whether $\texttt{COL_DISPERSION_RADIAL}$ is treated as component-, radial-, and/or section-independent. This field is optional. When left out, multiplexing behavior is inferred from the length of $\texttt{COL_DISPERSION_RADIAL}$. Valid modes are:
>
> 0. Component-independent, radial-independent, section-independent; length of $\texttt{COL_DISPERSION_RADIAL}$ is 1
> 1. Component-independent, radial-dependent, section-independent; length of $\texttt{COL_DISPERSION_RADIAL}$ is $\texttt{RAD_NELEM}$
> 2. Component-dependent, radial-independent, section-independent; length of $\texttt{COL_DISPERSION_RADIAL}$ is $\texttt{NCOMP}$
> 3. Component-dependent, radial-dependent, section-independent; length of $\texttt{COL_DISPERSION_RADIAL}$ is $\texttt{NCOMP} \cdot \texttt{RAD_NELEM}$; ordering is radial-major
> 4. Component-independent, radial-independent, section-dependent; length of $\texttt{COL_DISPERSION_RADIAL}$ is $\texttt{NSEC}$
> 5. Component-independent, radial-dependent, section-dependent; length of $\texttt{COL_DISPERSION_RADIAL}$ is $\texttt{RAD_NELEM} \cdot \texttt{NSEC}$; ordering is section-major
> 6. Component-dependent, radial-independent, section-dependent; length of $\texttt{COL_DISPERSION_RADIAL}$ is $\texttt{NCOMP} \cdot \texttt{NSEC}$; ordering is section-major
> 7. Component-dependent, radial-dependent, section-dependent; length of $\texttt{COL_DISPERSION_RADIAL}$ is $\texttt{NCOMP} \cdot \texttt{RAD_NELEM} \cdot \texttt{NSEC}$; ordering is section-radial-major
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\{0, \dots, 7 \}$
>   - **Length:** 1
> ```

`INIT_C`

> Initial concentrations for each component in all radial zones of the bulk mobile phase (length $\texttt{NCOMP}$), or for each component in each radial zone (length $\texttt{NCOMP} \cdot \texttt{RAD_NELEM}$, ordering radial-major)
>
> **Unit:** $\mathrm{mol}\,\mathrm{m}_{\mathrm{IV}}^{-3}$
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $\geq 0$
>   - **Length:** $\texttt{NCOMP} / \texttt{NCOMP} \cdot \texttt{RAD_NELEM}$
> ```

`INIT_STATE`

> Full state vector for initialization (optional, $\texttt{INIT_C}$, $\texttt{INIT_CP}$, and $\texttt{INIT_CS}$ will be ignored; if length is $2\texttt{NDOF}$, then the second half is used for time derivatives).
> The ordering of the state vector is defined in {ref}`unitoperationstateordering`.
>
> **Unit:** $various$
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $\mathbb{R}$
>   - **Length:** $\texttt{NDOF} / 2\texttt{NDOF}$
> ```

`NREAC_LIQUID`

> Number of liquid phase reaction models (optional, only if liquid reactions are present).
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq 0$
>   - **Length:** 1
> ```

## Group /input/model/unit_XXX/liquid_reaction_YYY

Each reaction is specified in another subgroup `liquid_reaction_YYY`, see {ref}`ffreaction`.

## Group /input/model/unit_XXX/particle_type_XXX

Each particle type is specified in another subgroup `particle_type_XXX`, see {ref}`particle-model-config`.

## Group /input/model/unit_XXX/discretization - UNIT_TYPE - COLUMN_MODEL_2D

`USE_ANALYTIC_JACOBIAN`

> Determines whether analytically computed Jacobian matrix (faster) is used (value is 1) instead of Jacobians generated by algorithmic differentiation (slower, value is 0)
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\{0, 1\}$
>   - **Length:** 1
> ```

## Spatial discretization - Numerical Methods

CADET offers two spatial discretization methods: Finite Volumes (FV) and Discontinuous Galerkin (DG). Each method has its own set of input fields.
While both methods approximate the same solution to the same underlying model, they may differ in terms of computational performance.
With our currently implemented variants of FV and DG, FV perform better for solutions with steep gradients or discontinuities, while DG can be much faster for rather smooth solutions.
For the same number of discrete points, DG will generally be slower but often more accurate.

For further information on the choice of discretization methods and their parameters, see {ref}`spatial-discretization-methods`.

`SPATIAL_METHOD`

> Spatial discretization method
>
> ```{list-table}
>
> * - **Type:** string
>   - **Range:** $\{\texttt{FV}, \texttt{DG}\}$
>   - **Length:** 1
> ```

`RADIAL_DISC_TYPE`

> Specifies the radial discretization scheme.
>
> ```{list-table}
>
> * - **Type:** string
>   - **Range:** $\{\texttt{EQUIDISTANT}, \texttt{EQUIVOLUME}, \texttt{USER_DEFINED} \}$
>   - **Length:** 1
> ```

`RADIAL_DISC_VECTOR`

> Coordinates for the radial compartment boundaries (ignored if $\texttt{RADIAL_DISC_TYPE} \neq \texttt{USER_DEFINED}$). The coordinates are absolute and have to include the endpoints 0 and $\texttt{COL_RADIUS}$. The values are expected in ascending order (i.e., 0 is the first and $\texttt{COL_RADIUS}$ the last value in the array).
>
> **Unit:** $\mathrm{m}$
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $[0,\texttt{COL_RADIUS}]$
>   - **Length:** $\texttt{RAD_NELEM} + 1$
> ```

## Discontinuous Galerkin

`AX_POLYDEG`

> DG polynomial degree for axial discretization. Optional, defaults to 4 and $N^z_d \in \{3, 4, 5\}$ is recommended.
> The total number of axial discrete points is given by (`AX_POLYDEG` + 1 ) * `AX_NELEM`
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq 1$
>   - **Length:** 1
> ```

`AX_NELEM`

> Number of axial column discretization DG elements. Optional and inferred from `NCOL` if left out. The total number of axial discrete points is given by (`AX_POLYDEG` + 1 ) * `AX_NELEM`
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq 1$
>   - **Length:** 1
> ```

`RAD_POLYDEG`

> DG polynomial degree for radial discretization. Optional, defaults to 4 and $N^r_d \in \{3, 4, 5\}$ is recommended, and should generally be the same as the axial degree.
> The total number of radial discrete points is given by (`RAD_POLYDEG` + 1 ) * `RAD_NELEM`
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq 1$
>   - **Length:** 1
> ```

`RAD_NELEM`

> Number of radial column discretization DG elements.
> Optional and inferred from `NRAD` if left out.
> The total number of radial discrete points is given by (`RAD_POLYDEG` + 1 ) * `RAD_NELEM`.
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq 1$
>   - **Length:** 1
> ```

`LINEAR_SOLVER`

> Specifies the linear solver variant used to factorize the semidiscretized system. Optional, defaults to `SparseLU`. For more information on these solvers, we refer to the [Eigen documentation](https://eigen.tuxfamily.org/)
>
> ```{list-table}
>
> * - **Type:** string
>   - **Range:** $\{\texttt{SparseLU}, \texttt{SparseQR}, ..., \texttt{BiCGSTAB}\}$
>   - **Length:** 1
> ```
>
> For further information on discretization parameters, see also {ref}`non-consistency-solver-parameters`.

## Finite Volumes

`NCOL`

> Number of axial column discretization cells
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq 1$
>   - **Length:** 1
> ```

`NRAD`

> Number of radial column discretization cells.
> Replaces the `RAD_NELEM` field in the length specification of parameters.
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq 1$
>   - **Length:** 1
> ```

`LINEAR_SOLVER_BULK`

> Linear solver used for the sparse column bulk block. This field is optional, the best available method is selected (i.e., sparse direct solver if possible). Valid values are:
>
> - $\texttt{DENSE}$ Converts the sparse matrix into a banded matrix and uses regular LAPACK. Slow and memory intensive, but always available.
> - $\texttt{UMFPACK}$ Uses the UMFPACK sparse direct solver (LU decomposition) from SuiteSparse. Fast, but has to be enabled when compiling and requires UMFPACK library.
> - $\texttt{SUPERLU}$ Uses the SuperLU sparse direct solver (LU decomposition). Fast, but has to be enabled when compiling and requires SuperLU library.
>
> ```{list-table}
>
> * - **Type:** string
>   - **Range:** $\{\texttt{DENSE},\texttt{UMFPACK},\texttt{SUPERLU}\}$
>   - **Length:** 1
> ```

`RECONSTRUCTION`

> Type of reconstruction method for fluxes
>
> ```{list-table}
>
> * - **Type:** string
>   - **Range:** $\texttt{WENO, KOREN}$
>   - **Length:** 1
> ```
>
> For further information on discretization parameters for reconstruction methods (including `AXIAL_GRID_FACES` for non-equidistant axial grid spacing), see also {ref}`flux-reconstruction-methods` (FV specific).

The following FV discretization parameters configure the GMRES solver of the unit operation's own
arrow-head Schur complement and are required for the FV discretization. They are independent of the
identically named fields in {ref}`ffmodelsolver`, which size the GMRES of the network.

`GS_TYPE`

> Type of Gram-Schmidt orthogonalization, see IDAS guide Section 4.5.7.3, p. 41f. A value of 0 enables classical Gram-Schmidt, a value of 1 uses modified Gram-Schmidt.
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\{0, 1\}$
>   - **Length:** 1
> ```

`MAX_KRYLOV`

> Defines the size of the Krylov subspace in the iterative linear GMRES solver (0: $\texttt{MAX_KRYLOV} = \texttt{NCOL} \cdot \texttt{NRAD} \cdot \texttt{NCOMP} \cdot \texttt{NPARTYPE}$)
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\{0, \dots, \texttt{NCOL} \cdot \texttt{NRAD} \cdot \texttt{NCOMP} \cdot \texttt{NPARTYPE} \}$
>   - **Length:** 1
> ```

`MAX_RESTARTS`

> Maximum number of restarts in the GMRES algorithm. If lack of memory is not an issue, better use a larger Krylov space than restarts.
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq 0$
>   - **Length:** 1
> ```

`SCHUR_SAFETY`

> Schur safety factor; Influences the tradeoff between linear iterations and nonlinear error control; see IDAS guide Section 2.1 and 5.
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $\geq 0$
>   - **Length:** 1
> ```

For further discretization parameters, see also {ref}`flux-reconstruction-methods`, and {ref}`non-consistency-solver-parameters`.
