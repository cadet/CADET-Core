(simplified-multi-state-steric-mass-action-config)=
# Simplified Multi-State Steric Mass Action

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = SIMPLIFIED_MULTISTATE_STERIC_MASS_ACTION**

For information on model equations, refer to {ref}`simplified-multi-state-steric-mass-action-model`.

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

`SMSSMA_LAMBDA`

: Stationary phase capacity (monovalent salt counterions); The total
  number of binding sites available on the resin surface

  **Unit:** $mol~m_{SP}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** 1
  ```

`SMSSMA_KA`

: Adsorption rate constants of the components to the different bound
  states in component-major ordering

  **Unit:** $m_{MP}^{3}~m_{SP}^{-3}~s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NTOTALBND
  ```

`SMSSMA_KD`

: Desorption rate constants of the components to the different bound
  states in component-major ordering

  **Unit:** $s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NTOTALBND
  ```

`SMSSMA_NU_MIN`

: Characteristic charges of the components in the first (weakest) bound
  state

  **Unit:** $s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NCOMP
  ```

`SMSSMA_NU_MAX`

: Characteristic charges of the components in the last (strongest)
  bound state

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NCOMP
  ```

`SMSSMA_NU_QUAD`

: Quadratic modifiers of the characteristic charges of the different
  components depending on the index of the bound state

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** NCOMP
  ```

`SMSSMA_SIGMA_MIN`

: Steric factors of the components in the first (weakest) bound state

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NCOMP
  ```

`SMSSMA_SIGMA_MAX`

: Steric factors of the components in the last (strongest) bound state

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NCOMP
  ```

`SMSSMA_SIGMA_QUAD`

: Quadratic modifiers of steric factors of the different components
  depending on the index of the bound state

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** NCOMP
  ```

`SMSSMA_KWS`

: Exchange rates from a weakly bound state to the next stronger bound
  state

  **Unit:** $s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NCOMP
  ```

`SMSSMA_KWS_LIN`

: Linear exchange rate coefficients from a weakly bound state to the
  next stronger bound state

  **Unit:** $s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** NCOMP
  ```

`SMSSMA_KWS_QUAD`

: Quadratic exchange rate coefficients from a weakly bound state to the
  next stronger bound state

  **Unit:** $s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** NCOMP
  ```

`SMSSMA_KSW`

: Exchange rates from a strongly bound state to the next weaker bound
  state

  **Unit:** $s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NCOMP
  ```

`SMSSMA_KSW_LIN`

: Linear exchange rate coefficients from a strongly bound state to the
  next weaker bound state

  **Unit:** $s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** NCOMP
  ```

`SMSSMA_KSW_QUAD`

: Quadratic exchange rate coefficients from a strongly bound state to
  the next weaker bound state

  **Unit:** $s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** NCOMP
  ```

`SMSSMA_NU_SALT`

: Characteristic charge of the salt component (optional, defaults to
  $1.0$)

  **Unit:** [-]

  ```{list-table}

  * - **Type:** double
    - **Range:** $\gt 0$
    - **Length:** 1
  ```

`SMSSMA_REFC0`

: Reference liquid phase concentration (optional, defaults to
  $1.0$)

  **Unit:** $mol~m_{MP}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\gt 0$
    - **Length:** 1
  ```

`SMSSMA_REFQ`

: Reference solid phase concentration (optional, defaults to
  $1.0$)

  **Unit:** $mol~m_{SP}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\gt 0$
    - **Length:** 1
  ```
