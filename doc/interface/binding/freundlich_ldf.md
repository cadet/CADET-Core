(freundlich-ldf-config)=
# Freundlich LDF

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = FREUNDLICH_LDF**

For information on model equations, refer to {ref}`freundlich-ldf-model`.

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

`FLDF_KKIN`

: Driving force coefficient for each component

  **Unit:** $s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** 1/NTOTALBND
  ```

`FLDF_KF`

: Freundlich coefficient for each component

  **Unit:** $m_{MP}^{3/n}~m_{SP}^{-3}~mol^{1-1/n}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** 1/NTOTALBND
  ```

`FLDF_N`

: Freundlich exponent for each component

  **Unit:** [-]

  ```{list-table}

  * - **Type:** double
    - **Range:** $> 0$
    - **Length:** 1/NTOTALBND
  ```
