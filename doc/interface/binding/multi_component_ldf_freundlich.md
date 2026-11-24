(multi-component-ldf-freundlich-config)=
# Multi Component Linear Driving Force Freundlich

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = MULTI_COMPONENT_LDF_FREUNDLICH**

For information on model equations, refer to {ref}`multi-component-ldf-freundlich-model`.

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

`MCLDFFRL_KLDF`

: Rate constants in linear driving force approach

  **Unit:** $s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NCOMP
  ```

`MCLDFFRL_KF`

: Proportionality constants

  **Unit:** $m_{MP}^{3/n}~m_{SP}^{-3}~mol^{1-1/n}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NCOMP
  ```

`MCLDFFRL_EXP`

: Freundlich exponent

  **Unit:** [-]

  ```{list-table}

  * - **Type:** double
    - **Range:** $\gt 0$
    - **Length:** NCOMP
  ```

`MCLDFFRL_A`

: Component influences in row-major ordering

  **Unit:** [-]

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** $\text{NCOMP}^2$
  ```

`MCLDFFRL_TAU`

: Small constant that ensures numerical stability

  **Unit:** $mol~m_{MP}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\gt 0$
    - **Length:** 1
  ```
