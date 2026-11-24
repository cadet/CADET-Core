(multi-component-sips-config)=
# Multi Component Sips

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = MULTI_COMPONENT_SIPS**

For information on model equations, refer to {ref}`multi-component-sips-model`.

`IS_KINETIC`

: Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
  quasi-stationary. If a single value is given, the mode is set for all
  bound states. Otherwise, the adsorption mode is set for each bound
  state separately.

  ```{list-table}

  * - **Type:** int
    - **Range:** {0,1}
    - **Length:** 1/NTOTALBND
  ```

`SIPS_KA`

: Adsorption rate constants

  **Unit:** $m_{MP}^3~mol^{-1}~s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NCOMP
  ```

`SIPS_KD`

: Desorption rate constants

  **Unit:** $s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NCOMP
  ```

`SIPS_QMAX`

: Maximum adsorption capacities

  **Unit:** $mol~m_{SP}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\gt 0$
    - **Length:** NCOMP
  ```

`SIPS_EXP`

: Freundlich-type exponent

  ```{list-table}

  * - **Type:** double
    - **Range:** $\gt 0$
    - **Length:** NCOMP
  ```

`SIPS_REFC0`

: Reference liquid phase concentration (optional, defaults to
  $1.0$)

  **Unit:** $mol~m_{MP}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\gt 0$
    - **Length:** 1
  ```

`SIPS_REFQ`

: Reference solid phase concentration (optional, defaults to
  $1.0$)

  **Unit:** $mol~m_{SP}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\gt 0$
    - **Length:** 1
  ```

`SIPS_LINEAR_THRESHOLD`

: Threshold for linearization. Originally, the liquid phase concentration
  enters the adsorption rate via a power term. Below this threshold, a
  quadratic approximation is used instead. This ensures numerical stability.

  **Unit:** $mol~m_{MP}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\gt 0$
    - **Length:** 1
  ```
