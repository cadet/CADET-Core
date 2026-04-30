(axial-flow-column-2d-config)=

# Axial Flow Column 2D

## Group /input/model/unit_XXX - UNIT_TYPE - COLUMN_MODEL_2D

`UNIT_TYPE`

> Specifies the type of unit operation model
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** string | **Range:** `COLUMN_MODEL_2D` | **Length:** 1 |

`NCOMP`

> Number of chemical components in the chromatographic medium
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 1$ | **Length:** 1 |

`CROSS_SECTION_AREA`

> Cross section area of the column. This parameter is optional and will be ignored if `COL_RADIUS` is provided.
> **Unit:** $\mathrm{m}^{2}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $>0$ | **Length:** 1 |

`COL_LENGTH`

> Column length / height
>
> **Unit:** $\mathrm{m}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $> 0$ | **Length:** 1 |

`COL_RADIUS`

> Column radius. This parameter is optional if `CROSS_SECTION_AREA` is provided.
>
> **Unit:** $\mathrm{m}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $> 0$ | **Length:** 1 |

`COL_POROSITY`

> Column porosity, either constant (length is 1) or for each radial zone (length is `NRAD`).
> In case of a spatially inhomogeneous setting, the `SENS_PARTYPE` field is used for indexing the radial zone when specifying parameter sensitivities.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $(0,1]$ | **Length:** $1 / \texttt{NRAD}$ |

`NPARTYPE`

> Number of particle types. Defaults to 0.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 0$ | **Length:** 1 |

`PAR_TYPE_VOLFRAC`

> Volume fractions of the particle types. The volume fractions can be set homogeneous or individually along both axes. For each cell, the volume fractions have to sum to 1.
> In case of a spatially inhomogeneous setting, the `SENS_SECTION` field is used for indexing the axial cell and the `SENS_REACTION` field is used for indexing the radial cell when specifying parameter sensitivities. This field is optional in case of only one particle type.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $[0,1]$ | **Length:** see `PAR_TYPE_VOLFRAC_MULTIPLEX` |

`PAR_TYPE_VOLFRAC_MULTIPLEX`

> > Multiplexing mode of `PAR_TYPE_VOLFRAC`. Determines whether `PAR_TYPE_VOLFRAC` is treated as radial- and/or section-independent.
> > This field is optional. When left out, multiplexing behavior is inferred from the length of `PAR_TYPE_VOLFRAC`. Valid modes are:
>
> 0. Radial-independent, axial-independent; length of `PAR_TYPE_VOLFRAC` is `NPARTYPE`
> 1. Radial-dependent, axial-independent; length of `PAR_TYPE_VOLFRAC` is $\texttt{NRAD} \cdot \texttt{NPARTYPE}$; ordering is radial-major
> 2. Axial-dependent; length of `PAR_TYPE_VOLFRAC` is $\texttt{NCOL} \cdot \texttt{NPARTYPE}$; ordering is axial-major
> 3. Radial-dependent, axial-dependent; length of `PAR_TYPE_VOLFRAC` is $\texttt{NCOL} \cdot \texttt{NRAD} \cdot \texttt{NPARTYPE}$; ordering is axial-radial-major

`VELOCITY`

> Indicates flow direction in each radial zone (forward if value is positive, backward if value is negative), see Section [](#MUOPGRMflow2D)). In case of a spatially inhomogeneous setting, the `SENS_PARTYPE` field is used for indexing the radial cell when specifying parameter sensitivities.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** see `VELOCITY_MULTIPLEX` |

`VELOCITY_MULTIPLEX`

> > Multiplexing mode of `VELOCITY`. Determines whether `VELOCITY` is treated as radial- and/or section-independent. This field is optional. When left out, multiplexing behavior is inferred from the length of `VELOCITY`. Valid modes are:
>
> 0. Radial-independent, section-independent; length of `VELOCITY` is 1
> 1. Radial-dependent, section-independent; length of `VELOCITY` is `NRAD`
> 2. Section-dependent; length of `VELOCITY` is `NSEC`
> 3. Radial-dependent, section-dependent; length of `VELOCITY` is $\texttt{NRAD} \cdot \texttt{NSEC}$; ordering is section-major

`COL_DISPERSION_AXIAL`

> Axial dispersion coefficient
>
> **Unit:** $\mathrm{m}_{\mathrm{IV}}^{2}\,\mathrm{s}^{-1}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** see `COL_DISPERSION_AXIAL_MULTIPLEX` |

`COL_DISPERSION_AXIAL_MULTIPLEX`

> > Multiplexing mode of `COL_DISPERSION_AXIAL`. Determines whether `COL_DISPERSION_AXIAL` is treated as component-, radial-, and/or section-independent. This field is optional. When left out, multiplexing behavior is inferred from the length of `COL_DISPERSION_AXIAL`. Valid modes are:
>
> 0. Component-independent, radial-independent, section-independent; length of `COL_DISPERSION_AXIAL` is 1
> 1. Component-independent, radial-dependent, section-independent; length of `COL_DISPERSION_AXIAL` is `NRAD`
> 2. Component-dependent, radial-independent, section-independent; length of `COL_DISPERSION_AXIAL` is `NCOMP`
> 3. Component-dependent, radial-dependent, section-independent; length of `COL_DISPERSION_AXIAL` is $\texttt{NCOMP} \cdot \texttt{NRAD}$; ordering is radial-major
> 4. Component-independent, radial-independent, section-dependent; length of `COL_DISPERSION_AXIAL` is `NSEC`
> 5. Component-independent, radial-dependent, section-dependent; length of `COL_DISPERSION_AXIAL` is $\texttt{NRAD} \cdot \texttt{NSEC}$; ordering is section-major
> 6. Component-dependent, radial-independent, section-independent; length of `COL_DISPERSION_AXIAL` is $\texttt{NCOMP} \cdot \texttt{NSEC}$; ordering is section-major
> 7. Component-dependent, radial-dependent, section-dependent; length of `COL_DISPERSION_AXIAL` is $\texttt{NCOMP} \cdot \texttt{NRAD} \cdot \texttt{NSEC}$; ordering is section-radial-major
>
> > |   |   |   |
> > | --- | --- | --- |
> > | **Type:** int | **Range:** $\{0, \dots, 7 \}$ | **Length:** 1 |

`COL_DISPERSION_RADIAL`

> Radial dispersion coefficient. In case of a spatially inhomogeneous setting, the `SENS_PARTYPE` field is used for indexing the radial zone when specifying parameter sensitivities.
>
> **Unit:** $\mathrm{m}_{\mathrm{IV}}^{2}\,\mathrm{s}^{-1}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** see `COL_DISPERSION_RADIAL_MULTIPLEX` |

`COL_DISPERSION_RADIAL_MULTIPLEX`

> > Multiplexing mode of `COL_DISPERSION_RADIAL`. Determines whether `COL_DISPERSION_RADIAL` is treated as component-, radial-, and/or section-independent. This field is optional. When left out, multiplexing behavior is inferred from the length of `COL_DISPERSION_RADIAL`. Valid modes are:
>
> 0. Component-independent, radial-independent, section-independent; length of `COL_DISPERSION_RADIAL` is 1
> 1. Component-independent, radial-dependent, section-independent; length of `COL_DISPERSION_RADIAL` is `NRAD`
> 2. Component-dependent, radial-independent, section-independent; length of `COL_DISPERSION_RADIAL` is `NCOMP`
> 3. Component-dependent, radial-dependent, section-independent; length of `COL_DISPERSION_RADIAL` is $\texttt{NCOMP} \cdot \texttt{NRAD}$; ordering is radial-major
> 4. Component-independent, radial-independent, section-dependent; length of `COL_DISPERSION_RADIAL` is `NSEC`
> 5. Component-independent, radial-dependent, section-dependent; length of `COL_DISPERSION_RADIAL` is $\texttt{NRAD} \cdot \texttt{NSEC}$; ordering is section-major
> 6. Component-dependent, radial-independent, section-independent; length of `COL_DISPERSION_RADIAL` is $\texttt{NCOMP} \cdot \texttt{NSEC}$; ordering is section-major
> 7. Component-dependent, radial-dependent, section-dependent; length of `COL_DISPERSION_RADIAL` is $\texttt{NCOMP} \cdot \texttt{NRAD} \cdot \texttt{NSEC}$; ordering is section-radial-major
>
> > |   |   |   |
> > | --- | --- | --- |
> > | **Type:** int | **Range:** $\{0, \dots, 7 \}$ | **Length:** 1 |

`INIT_C`

> Initial concentrations for each component in all radial zones the bulk mobile phase (length `NCOMP`), or for each component in each radial zone (length $\texttt{NCOMP} \cdot \texttt{NRAD}$, ordering radial-major)
>
> **Unit:** $\mathrm{mol}\,\mathrm{m}_{\mathrm{IV}}^{-3}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ |  |

`INIT_STATE`

> Full state vector for initialization (optional, `INIT_C`, `INIT_CP`, and `INIT_CS` will be ignored; if length is $2\texttt{NDOF}$, then the second half is used for time derivatives).
> The ordering of the state vector is defined in [](#UnitOperationStateOrdering).
>
> **Unit:** $various$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** $\texttt{NDOF} / 2\texttt{NDOF}$ |

`NREAC_LIQUID`

> Number of liquid phase reaction models (optional, only if liquid reactions are present).
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 0$ | **Length:** 1 |

## Group /input/model/unit_XXX/liquid_reaction_YYY

Each reaction is specified in another subgroup `liquid_reaction_YYY`, see [](#FFReaction).

## Group /input/model/unit_XXX/particle_type_XXX

Each particle type is specified in another subgroup `particle_type_XXX`, see [](#particle-model-config).

## Group /input/model/unit_XXX/discretization - UNIT_TYPE - COLUMN_MODEL_2D

`USE_ANALYTIC_JACOBIAN`

> Determines whether analytically computed Jacobian matrix (faster) is used (value is 1) instead of Jacobians generated by algorithmic differentiation (slower, value is 0)
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\{0, 1\}$ | **Length:** 1 |

## Spatial discretization - Numerical Methods

CADET offers two spatial discretization methods: Finite Volumes (FV) and Discontinuous Galerkin (DG). Each method has it's own set of input fields.
While both methods approximate the same solution to the same underlying model, they may differ in terms of computational performance.
With our currently implemented variants of FV and DG, FV perform better for solutions with steep gradients or discontinuities, while DG can be much faster for rather smooth solutions.
For the same number of discrete points, DG will generally be slower but often more accurate.

For further information on the choice of discretization methods and their parameters, see [](#spatial-discretization-methods).

`SPATIAL_METHOD`

> Spatial discretization method. Optional, defaults to `FV`
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** string | **Range:** $\{\texttt{FV}, \texttt{DG}\}$ | **Length:** 1 |

`RADIAL_DISC_TYPE`

> Specifies the radial discretization scheme. Valid values are `EQUIDISTANT`, `EQUIVOLUME`, and `USER_DEFINED`.
>
> |   |   |
> | --- | --- |
> | **Type:** string | **Length:** 1 |

`RADIAL_DISC_VECTOR`

> Coordinates for the radial compartment boundaries (ignored if $\texttt{RADIAL\_DISC\_TYPE} \neq \texttt{USER\_DEFINED}$). The coordinates are absolute and have to include the endpoints 0 and `COLUMN_RADIUS`. The values are expected in ascending order (i.e., 0 is the first and `COLUMN_RADIUS` the last value in the array).
>
> **Unit:** $\mathrm{m}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $[0,\texttt{COLUMN\_RADIUS}]$ | **Length:** :math:`\texttt{NRAD} + 1 |

## Discontinuous Galerkin

`AX_POLYDEG`

> DG polynomial degree for axial discretization. Optional, defaults to 4 and $N^z_d \in \{3, 4, 5\}$ is recommended.
> The total number of axial discrete points is given by (`AX_POLYDEG` + 1 ) * `AX_NELEM`
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 1$ | **Length:** 1 |

`AX_NELEM`

> Number of axial column discretization DG cellselements. The total number of axial discrete points is given by (`POLYDEG` + 1 ) * `NELEM`
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 1$ | **Length:** 1 |

`RAD_POLYDEG`

> DG polynomial degree for radial discretization. Optional, defaults to 4 and $N^r_d \in \{3, 4, 5\}$ is recommended, and should generally be the same as the axial degree.
> The total number of radial discrete points is given by (`RAD_POLYDEG` + 1 ) * `RAD_NELEM`
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 1$ | **Length:** 1 |

`RAD_NELEM`

> Number of radial column discretization DG cellselements. The total number of axial discrete points is given by (`POLYDEG` + 1 ) * `NELEM`
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 1$ | **Length:** 1 |

`LINEAR_SOLVER`

> Specifies the linear solver variant used to factorize the semidiscretized system. Optional, defaults to `SparseLU`. For more information on these solvers, we refer to the [Eigen documentation](https://eigen.tuxfamily.org/)
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\{\texttt{SparseLU}, \texttt{SparseQR}, ..., \texttt{BiCGSTAB}\}$ | **Length:** 1 |
>
> For further information on discretization parameters, see also [](#non-consistency-solver-parameters).

## Finite Volumes

`NCOL`

> Number of axial column discretization cells
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 1$ | **Length:** 1 |

`NRAD`

> Number of radial column discretization cells
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 1$ | **Length:** 1 |

`NPAR`

> Number of particle (radial) discretization cells for each particle type
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 0$ | **Length:** $1 / \texttt{NPARTYPE}$ |

`LINEAR_SOLVER_BULK`

> > Linear solver used for the sparse column bulk block. This field is optional, the best available method is selected (i.e., sparse direct solver if possible). Valid values are:
>
> - `DENSE` Converts the sparse matrix into a banded matrix and uses regular LAPACK. Slow and memory intensive, but always available.
> - `UMFPACK` Uses the UMFPACK sparse direct solver (LU decomposition) from SuiteSparse. Fast, but has to be enabled when compiling and requires UMFPACK library.
> - `SUPERLU` Uses the SuperLU sparse direct solver (LU decomposition). Fast, but has to be enabled when compiling and requires SuperLU library.
>
> > |   |   |   |
> > | --- | --- | --- |
> > | **Type:** string | **Range:** $\{\texttt{DENSE},\texttt{UMFPACK},\texttt{SUPERLU}\}$ | **Length:** 1 |

`RECONSTRUCTION`

> Type of reconstruction method for fluxes
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** string | **Range:** `WENO, KOREN` | **Length:** 1 |
>
> For further information on discretization parameters for reconstruction methods (including `AXIAL_GRID_FACES` for non-equidistant axial grid spacing), see also [](#flux-reconstruction-methods) (FV specific).

`GS_TYPE`

> Type of Gram-Schmidt orthogonalization, see IDAS guide Section~4.5.7.3, p.~41f. A value of 0 enables classical Gram-Schmidt, a value of 1 uses modified Gram-Schmidt.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\{0, 1\}$ | **Length:** 1 |

`MAX_KRYLOV`

> Defines the size of the Krylov subspace in the iterative linear GMRES solver (0: $\texttt{MAX\_KRYLOV} = \texttt{NCOL} \cdot \texttt{NRAD} \cdot \texttt{NCOMP} \cdot \texttt{NPARTYPE}$)
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\{0, \dots, \texttt{NCOL} \cdot \texttt{NRAD} \cdot \texttt{NCOMP} \cdot \texttt{NPARTYPE} \}$ | **Length:** 1 |

`MAX_RESTARTS`

> Maximum number of restarts in the GMRES algorithm. If lack of memory is not an issue, better use a larger Krylov space than restarts.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 0$ | **Length:** 1 |

`SCHUR_SAFETY`

> Schur safety factor; Influences the tradeoff between linear iterations and nonlinear error control; see IDAS guide Section~2.1 and 5.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** 1 |

For further discretization parameters, see also [](#flux-reconstruction-methods), and [](#non-consistency-solver-parameters).
