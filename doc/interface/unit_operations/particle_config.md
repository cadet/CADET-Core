(particle-model-config)=

# Particle Model

## Group /input/model/unit_XXX/particle_type_XXX

`PAR_POROSITY`

> Particle porosity of all particle types or for each particle type
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $(0,1]$
>   - **Length:** 1
> ```

`PAR_RADIUS`

> Particle radius of all particle types or for each particle type
>
> **Unit:** $\mathrm{m}$
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $>0$
>   - **Length:** 1
> ```

`PAR_CORERADIUS`

> Particle core radius of all particle types or for each particle type (optional, defaults to $\mathrm{0}$)
> Is only applied when $\texttt{HAS_PORE_DIFFUSION} == 1$.
>
> **Unit:** $\mathrm{m}$
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $[0, \texttt{PAR_RADIUS})$
>   - **Length:** 1
> ```

`PORE_ACCESSIBILITY`

> Pore accessibility factor of each component in each particle type (optional, defaults to $1$).
> Note: When combined with a binding model, components that compete for the same bound states must use the same value, see {ref}`muopgrmsizeexclusion`.
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $(0, 1]$
>   - **Length:** see $\texttt{PORE_ACCESSIBILITY_MULTIPLEX}$
> ```

`PORE_ACCESSIBILITY_MULTIPLEX`

> Multiplexing mode of $\texttt{PORE_ACCESSIBILITY}$. Determines whether $\texttt{PORE_ACCESSIBILITY}$ is treated as component-, type-, and/or section-independent. This field is optional. When left out, multiplexing behavior is inferred from the length of $\texttt{PORE_ACCESSIBILITY}$. Valid modes are:
>
> 0. Component-dependent, section-independent; length of $\texttt{PORE_ACCESSIBILITY}$ is $\texttt{NCOMP}$
> 1. Component-dependent, section-dependent; length of $\texttt{PORE_ACCESSIBILITY}$ is $\texttt{NCOMP} \cdot \texttt{NSEC}$; ordering is section-major
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\{0, 1 \}$
>   - **Length:** 1
> ```

`HAS_FILM_DIFFUSION`

> Specifies whether transport into the particles is limited by film diffusion kinetics.
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\{0, 1 \}$
>   - **Length:** 1
> ```

`FILM_DIFFUSION`

> Film diffusion coefficients for each component of each particle type, required if $\texttt{HAS_FILM_DIFFUSION} == 1$, otherwise ignored.
>
> **Unit:** $\mathrm{m}\,\mathrm{s}^{-1}$
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $\geq 0$
>   - **Length:** see $\texttt{FILM_DIFFUSION_MULTIPLEX}$
> ```

`FILM_DIFFUSION_MULTIPLEX`

> Multiplexing mode of $\texttt{FILM_DIFFUSION}$. Determines whether $\texttt{FILM_DIFFUSION}$ is treated as component-, type-, and/or section-independent. This field is optional. When left out, multiplexing behavior is inferred from the length of $\texttt{FILM_DIFFUSION}$. Valid modes are:
>
> 0. Component-dependent, section-independent; length of $\texttt{FILM_DIFFUSION}$ is $\texttt{NCOMP}$
> 1. Component-dependent, section-dependent; length of $\texttt{FILM_DIFFUSION}$ is $\texttt{NCOMP} \cdot \texttt{NSEC}$; ordering is section-major
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\{0, 1 \}$
>   - **Length:** 1
> ```

`HAS_PORE_DIFFUSION`

> Specifies whether radial transport within the particle pores is limited by pore diffusion kinetics. Optional, defaults to $0$.
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\{0, 1 \}$
>   - **Length:** 1
> ```

`PORE_DIFFUSION`

> Effective particle diffusion coefficients of each component in each particle type.
> The field is still read when $\texttt{HAS_PORE_DIFFUSION} == 0$, where its value has no
> effect on the model and is conventionally given as $0$.
>
> **Unit:** $\mathrm{m}_{\mathrm{MP}}^{2}\,\mathrm{s}^{-1}$
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $\geq 0$
>   - **Length:** see $\texttt{PORE_DIFFUSION_MULTIPLEX}$
> ```

`PORE_DIFFUSION_MULTIPLEX`

> Multiplexing mode of $\texttt{PORE_DIFFUSION}$. Determines whether $\texttt{PORE_DIFFUSION}$ is treated as component-, type-, and/or section-independent. This field is optional. When left out, multiplexing behavior is inferred from the length of $\texttt{PORE_DIFFUSION}$. Valid modes are:
>
> 0. Component-dependent, section-independent; length of $\texttt{PORE_DIFFUSION}$ is $\texttt{NCOMP}$
> 1. Component-dependent, section-dependent; length of $\texttt{PORE_DIFFUSION}$ is $\texttt{NCOMP} \cdot \texttt{NSEC}$; ordering is section-major
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\{0, 1 \}$
>   - **Length:** 1
> ```

`HAS_SURFACE_DIFFUSION`

> Specifies whether radial transport within the particle is supported but limited by surface diffusion kinetics. Optional, defaults to $0$.
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\{0, 1 \}$
>   - **Length:** 1
> ```

`SURFACE_DIFFUSION`

> Particle surface diffusion coefficients of each bound state of each component in each particle type.
> Optional, defaults to all 0 $\mathrm{m}_{\mathrm{SP}}^{2}\,\mathrm{s}^{-1}$, required if $\texttt{HAS_SURFACE_DIFFUSION} == 1$, otherwise ignored.
>
> **Unit:** $\mathrm{m}_{\mathrm{SP}}^{2}\,\mathrm{s}^{-1}$
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $\geq 0$
>   - **Length:** see $\texttt{SURFACE_DIFFUSION_MULTIPLEX}$
> ```

`SURFACE_DIFFUSION_MULTIPLEX`

: Multiplexing mode of $\texttt{SURFACE_DIFFUSION}$. Determines whether $\texttt{SURFACE_DIFFUSION}$ is treated as component-, type-, and/or section-independent. This field is optional. When left out, multiplexing behavior is inferred from the length of $\texttt{SURFACE_DIFFUSION}$. Valid modes are:

  0. Component-dependent, section-independent; length of $\texttt{SURFACE_DIFFUSION}$ is $\texttt{NBOUND}$
  1. Component-dependent, section-dependent; length of $\texttt{SURFACE_DIFFUSION}$ is $\texttt{NBOUND} \cdot \texttt{NSEC}$; ordering is section-major

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0, 1 \}$
    - **Length:** 1
  ```

`PAR_GEOM`

> Specifies the particle geometry for all or each particle type.
> Valid values are $\texttt{SPHERE}$, $\texttt{CYLINDER}$, $\texttt{SLAB}$. Optional, defaults to $\texttt{SPHERE}$.
> Is only applied when $\texttt{HAS_PORE_DIFFUSION} == 1$.
>
> ```{list-table}
>
> * - **Type:** string
>   - **Range:** $\{\texttt{SPHERE}, \texttt{CYLINDER}, \texttt{SLAB} \}$
>   - **Length:** 1
> ```

`ADSORPTION_MODEL`

> Specifies the type of binding model of each particle type (or of all particle types if length is $1$)
>
> ```{list-table}
>
> * - **Type:** string
>   - **Range:** {ref}`ffadsorption`
>   - **Length:** 1
> ```

`NBOUND`

> Number of bound states for each component of this particle type. The bound states are kept if $\texttt{ADSORPTION_MODEL}$ is $\texttt{NONE}$, which leaves states that solid phase and cross-phase reactions can act on.
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq 0$
>   - **Length:** $\texttt{NCOMP}$
> ```

`INIT_CP`

> Initial concentrations for each component in the bead liquid phase (optional, $\texttt{INIT_C}$ is used if left out).
> The length of this field is $\texttt{NCOMP}$.
> Only in case of a 2D bulk model, the field length *can* also be $\texttt{NCOMP} \cdot \texttt{NRAD}$ to specify radial dependence, in which case the ordering is radial position major.
> If $\texttt{BINDING_PARTYPE_DEPENDENT}$ is $0$, the values across different particle types (if multiple particle types are being used) must be the same.
>
> **Unit:** $\mathrm{mol}\,\mathrm{m}_{\mathrm{MP}}^{-3}$
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $\geq 0$
>   - **Length:** $\texttt{NCOMP}$
> ```

`INIT_CS`

> Initial concentrations for each bound state of each component in the bead solid phase.
> The length of this field is $\texttt{NBOUND}$.
> Only in case of a 2D bulk model, the field length *can* also be $\texttt{NBOUND} \cdot \texttt{NRAD}$ to specify radial dependence, in which case the ordering is radial position major.
> If $\texttt{BINDING_PARTYPE_DEPENDENT}$ is $0$, the values across different particle types (if multiple particle types are being used) must be the same.
>
> **Unit:** $\mathrm{mol}\,\mathrm{m}_{\mathrm{SP}}^{-3}$
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $\geq 0$
>   - **Length:** $\texttt{NBOUND}$
> ```

`PARAMNAME_PARTYPE_DEPENDENT`

> Only required for parameter sensitivities of the respective parameter, defaults to 1.
> Specifies whether or not a parameter (parameter name substitutes the `PARAMNAME` prefix of field name) is the same across particle types or 'dependent' on the particle type.
> In the first case, the parameter sensitivity should be computed jointly across particle types, by specifying this field as 0.
> Can be specified for any of the above parameters except `NCOMP` and `NBOUND`.
> For more information on parameter sensitivities, see {ref}`sensitivity`.
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\{0, 1 \}$
>   - **Length:** 1
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

`NREAC_SOLID`

: Number of solid phase reaction models (optional, only if solid reactions are present).

  ```{list-table}

  * - **Type:** int
    - **Range:** $\geq 0$
    - **Length:** 1
  ```

`NREAC_CROSS_PHASE`

: Number of cross-phase reaction models (optional, only if cross-phase reactions are present).

  ```{list-table}

  * - **Type:** int
    - **Range:** $\geq 0$
    - **Length:** 1
  ```

## Group /input/model/unit_XXX/phase_reaction_YYY

Each reaction is specified in another subgroup `phase_reaction_YYY`, see {ref}`ffreaction`.

## Group /input/model/unit_XXX/particle_type_XXX/discretization

`PAR_DISC_TYPE`

> Specifies the discretization scheme inside the particles for all or each particle type. Optional, defaults to $\texttt{EQUIDISTANT}$.
> Any other value is rejected with an error.
>
> ```{list-table}
>
> * - **Type:** string
>   - **Range:** $\{\texttt{EQUIDISTANT}, \texttt{EQUIVOLUME}, \texttt{USER_DEFINED} \}$
>   - **Length:** 1
> ```

`PAR_DISC_VECTOR`

> Node coordinates for the DG element (or FV cell) boundaries (ignored if $\texttt{PAR_DISC_TYPE} \neq \texttt{USER_DEFINED}$).
> The coordinates are relative and have to include the endpoints $0$ and $1$.
> They are later linearly mapped to the true radial range $[r_{c,j}, r_{p,j}]$.
> The coordinates for each particle type are appended to one long vector in particle type major ordering.
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $[0,1]$
>   - **Length:** $\texttt{PAR_NELEM} + 1$ (DG) / $\texttt{NCELLS} + 1$ (FV)
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

## Discontinuous Galerkin

`PAR_POLYDEG`

> DG particle (radial) polynomial degree. The total number of particle (radial) discrete points is given by (`PAR_POLYDEG` + 1 ) * `PAR_NELEM`.
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq 1$
>   - **Length:** 1
> ```

`PAR_NELEM`

> Number of particle (radial) discretization DG elements for each particle type. For the particle discretization, it is usually most performant to fix `PAR_NELEM` = 1 and to increase the polynomial degree for more accuracy.
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq 1$
>   - **Length:** 1
> ```

`PAR_GSM`

> Specifies whether Galerkin spectral method should be used (as opposed to discontinuous variant, DGSEM), optional, defaults to 1 if `PAR_NELEM` == 1. Always recommended.
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\{ 0,1 \}$
>   - **Length:** 1
> ```
>
> For further discretization parameters, see also {ref}`non-consistency-solver-parameters`.

## Finite Volumes

`NCELLS`

> Number of particle (radial) discretization points for each particle type
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq 1$
>   - **Length:** 1
> ```

`FV_BOUNDARY_ORDER`

> Order of accuracy of outer particle boundary condition. Optional, defaults to $2$.
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\{ 1,2 \}$
>   - **Length:** 1
> ```

`OPTIMIZE_PAR_BANDWIDTH`

> Only available if both particle and bulk spatial methods are FV, and specified in the unit operation's discretization group, i.e. in `/input/model/unit_XXX/discretization`.
> Determines whether the surface diffusion parameters $\texttt{SURFACE_DIFFUSION}$ are fixed if the parameters are zero.
> If the parameters are fixed to zero ($\texttt{OPTIMIZE_PAR_BANDWIDTH} = 1$, $\texttt{SURFACE_DIFFUSION} = 0$), the parameters must not become non-zero during this or subsequent simulation runs.
> The internal data structures are optimized for a more efficient simulation. This field is optional and defaults to $1$ (optimization enabled).
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\{0, 1\}$
>   - **Length:** 1
> ```

When using the FV method, we generally recommend specifying `USE_MODIFIED_NEWTON = 0` in {ref}`ffsolvertime`, i.e. to use the full Newton method to solve the linear system within the time integrator.
For further discretization parameters, see also {ref}`flux-reconstruction-methods` (FV specific), and {ref}`non-consistency-solver-parameters`.
