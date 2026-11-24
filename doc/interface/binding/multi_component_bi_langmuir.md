(multi-component-bi-langmuir-config)=
# Multi Component Bi-Langmuir

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = MULTI_COMPONENT_BILANGMUIR**

For information on model equations, refer to {ref}`multi-component-bi-langmuir-model`.

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

`MCBL_KA`

: Adsorption rate constants in state-major ordering (see {ref}`ordering-multi-dimensional-data`)

  **Unit:** $m_{MP}^3~mol^{-1}~s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NSTATES $\cdot$ NCOMP
  ```

`MCBL_KD`

: Desorption rate constants in state-major ordering

  **Unit:** $s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NSTATES $\cdot$ NCOMP
  ```

`MCBL_QMAX`

: Maximum adsorption capacities in state-major ordering

  **Unit:** $mol~m_{SP}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NSTATES $\cdot$ NCOMP
  ```
