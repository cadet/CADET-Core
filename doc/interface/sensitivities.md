(sensitivity)=

# Parameter Sensitivities

(ffsensitivity)=

## Group /input/sensitivity

`NSENS`

> Number of sensitivities to be computed
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq 0$
>   - **Length:** 1
> ```

`SENS_METHOD`

> Method used for computation of sensitivities (algorithmic differentiation)
>
> ```{list-table}
>
> * - **Type:** string
>   - **Range:** $\texttt{ad1}$
>   - **Length:** 1
> ```

(ffsensitivityparam)=

## Group /input/sensitivity/param_XXX

`SENS_UNIT`

> Unit operation index
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq 0$
>   - **Length:** $\geq 1$
> ```

`SENS_NAME`

> Name of the parameter. Some parameters do not support sensitivities with every spatial discretization, see {ref}`unsupported-sensitivities` below. Requesting one of those is rejected with an error rather than silently answered with a wrong number.
>
> ```{list-table}
>
> * - **Type:** string
>   - **Length:** $\geq 1$
> ```

`SENS_COMP`

> Component index ($-1$ if parameter is independent of components, see the multiplexing of the corresponding parameter, e.g. `FILM_DIFFUSION_MULTIPLEX`)
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq -1$
>   - **Length:** $\geq 1$
> ```

`SENS_PARTYPE`

> Particle type index ($-1$ if parameter is independent of particle types, see the multiplexing of the corresponding parameter, e.g. `PORE_DIFFUSION_MULTIPLEX`)
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq -1$
>   - **Length:** $\geq 1$
> ```

`SENS_REACTION`

> Reaction index ($-1$ if parameter is independent of reactions)
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq -1$
>   - **Length:** $\geq 1$
> ```

`SENS_BOUNDPHASE`

> Bound phase index ($-1$ if parameter is independent of bound phases)
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq -1$
>   - **Length:** $\geq 1$
> ```

`SENS_SECTION`

> Section index ($-1$ if parameter is independent of sections)
>
> ```{list-table}
>
> * - **Type:** int
>   - **Range:** $\geq -1$
>   - **Length:** $\geq 1$
> ```

`SENS_ABSTOL`

> Absolute tolerance used in the computation of the sensitivities (optional, defaults to $10^{-5}$). Rule of thumb: $\texttt{ABSTOL} / \texttt{PARAM_VALUE}$
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $\geq 0.0$
>   - **Length:** $\geq 1$
> ```

`SENS_FACTOR`

> Linear factor of the combined sensitivity (optional, taken as $1.0$ if left out)
>
> ```{list-table}
>
> * - **Type:** double
>   - **Range:** $\mathbb{R}$
>   - **Length:** $\geq 1$
> ```

(unsupported-sensitivities)=

## Unsupported sensitivities

A parameter enters a spatial operator either through the residual, where the algorithmic
differentiation carries its derivative along, or through an operator matrix that is assembled once
in double precision, where the derivative is lost. The second case cannot be answered correctly and
is rejected with an error. The alternatives below are the ones that have been checked against
central finite differences.

```{list-table}
:header-rows: 1

* - Parameter
  - Not available for
  - Alternative
* - `SECTION_TIMES`
  - every unit operation and discretization
  - Approximate with finite differences, i.e. rerun the simulation with a shifted section time. The derivative is dominated by the jump the solution makes when the switch moves, which needs a sensitivity jump condition that is not implemented.
* - `PAR_RADIUS`, `PAR_CORERADIUS`
  - DG discretized particles
  - Use the finite volume particle discretization.
* - `COL_DISPERSION`, `BED_LENGTH`, `COL_RADIUS_SMALL_END`, `COL_RADIUS_LARGE_END`
  - exact integration DG (`USE_COLLOCATION_DG = 0`) of an axial flow cylinder
  - Use collocation DG (`USE_COLLOCATION_DG = 1`) or the finite volume discretization. Those describe the cylinder through `CROSS_SECTION_AREA` rather than the two end radii.
* - `COL_DISPERSION`
  - DG discretized radial flow and frustum columns
  - Use the finite volume discretization.
* - `BED_LENGTH`
  - DG discretized radial flow and frustum columns
  - Use the finite volume discretization, which rebuilds its cell geometry when the bed length changes. This requires an equidistant grid; a `GRID_FACES` grid prescribes the cell faces independently of the bed length and rejects this sensitivity as well.
* - `COL_RADIUS_SMALL_END`, `COL_RADIUS_LARGE_END`
  - DG discretized radial flow and frustum columns
  - Approximate with finite differences. The finite volume discretization of these geometries derives the end radii from the cross section areas and has no such parameter.
* - `COL_DISPERSION_RADIAL`, `COL_RADIUS`, `BED_LENGTH`
  - DG discretized 2D column model
  - Approximate with finite differences. Note that `COL_DISPERSION_AXIAL` *is* available for this model.
```

See the issues [#214](https://github.com/cadet/CADET-Core/issues/214) and
[#810](https://github.com/cadet/CADET-Core/issues/810) for the background.
