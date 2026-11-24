(mobile-phase-modulator-langmuir-config)=
# Mobile Phase Modulator Langmuir

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = MOBILE_PHASE_MODULATOR**

For information on model equations, refer to {ref}`mobile-phase-modulator-langmuir-model`.

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

`MPM_KA`

: Adsorption rate constants

  **Unit:** $m_{MP}^3~mol^{-1}~s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NCOMP
  ```

`MPM_KD`

: Desorption rate constants

  **Unit:** $m_{MP}^{3\beta}~mol^{-\beta}~s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NCOMP
  ```

`MPM_QMAX`

: Maximum adsorption capacities

  **Unit:** $mol~m_{SP}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NCOMP
  ```

`MPM_BETA`

: Parameters describing the ion-exchange characteristics (IEX)

  **Unit:** $mol~m_{SP}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NCOMP
  ```

`MPM_GAMMA`

: Parameters describing the hydrophobicity (HIC)

  **Unit:** $m_{MP}^{3} mol^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** NCOMP
  ```

`MPM_LINEAR_THRESHOLD`

: Concentration of the modifier component below which the exponential
  term is linearized (optional)

  **Unit:** $mol~m_{MP}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** 1
  ```
