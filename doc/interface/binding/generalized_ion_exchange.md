(generalized-ion-exchange-config)=
# Generalized Ion Exchange

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = GENERALIZED_ION_EXCHANGE**

For information on model equations, refer to {ref}`generalized-ion-exchange-model`.

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

`GIEX_KA`

: Base value of adsorption rate constant

  **Unit:** $m_{MP}^{3}~m_{SP}^{-3}~s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NCOMP
  ```

`GIEX_KA_LIN`

: Coefficient of linear dependence of adsorption rate constant on
  modifier component

  **Unit:** $\text{[Mod]}^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Length:** NCOMP
  ```

`GIEX_KA_QUAD`

: Coefficient of quadratic dependence of adsorption rate constant on
  modifier component

  **Unit:** $\text{[Mod]}^{-2}$

  ```{list-table}

  * - **Type:** double
    - **Length:** NCOMP
  ```

`GIEX_KA_SALT`

: Salt coefficient of adsorption rate constants; difference of
  water-protein and salt-protein interactions

  ```{list-table}

  * - **Type:** double
    - **Length:** NCOMP
  ```

`GIEX_KA_PROT`

: Protein coefficient of adsorption rate constants; difference of
  water-protein and protein-protein interactions

  **Unit:** $m_{MP}^{3} mol^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Length:** NCOMP
  ```

`GIEX_KD`

: Base value of desorption rate constant

  **Unit:** $s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NCOMP
  ```

`GIEX_KD_LIN`

: Coefficient of linear dependence of desorption rate constant on
  modifier component

  **Unit:** $\text{[Mod]}^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Length:** NCOMP
  ```

`GIEX_KD_QUAD`

: Coefficient of quadratic dependence of desorption rate constant on
  modifier component

  **Unit:** $\text{[Mod]}^{-2}$

  ```{list-table}

  * - **Type:** double
    - **Length:** NCOMP
  ```

`GIEX_KD_SALT`

: Salt coefficient of desorption rate constants; difference of
  water-protein and salt-protein interactions

  ```{list-table}

  * - **Type:** double
    - **Length:** NCOMP
  ```

`GIEX_KD_PROT`

: Protein coefficient of desorption rate constants; difference of
  water-protein and protein-protein interactions

  **Unit:** $m_{MP}^{3} mol^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Length:** NCOMP
  ```

`GIEX_NU`

: Base value for characteristic charges of the protein; The number of
  sites $\nu$ that the protein interacts with on the resin
  surface

  ```{list-table}

  * - **Type:** double
    - **Length:** NCOMP
  ```

`GIEX_NU_LIN`

: Coefficient of linear dependence of characteristic charge on modifier
  component

  **Unit:** $\text{[Mod]}^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Length:** NCOMP
  ```

`GIEX_NU_QUAD`

: Coefficient of quadratic dependence of characteristic charge on
  modifier component

  **Unit:** $\text{[Mod]}^{-2}$

  ```{list-table}

  * - **Type:** double
    - **Length:** NCOMP
  ```

`GIEX_SIGMA`

: Steric factors of the protein; The number of sites $\sigma$ on
  the surface that are shielded by the protein and prevented from
  exchange with the salt counterions in solution

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** NCOMP
  ```

`GIEX_LAMBDA`

: Stationary phase capacity (monovalent salt counterions); The total
  number of binding sites available on the resin surface

  **Unit:** $mol~m_{SP}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** 1
  ```

`GIEX_REFC0`

: Reference liquid phase concentration (optional, defaults to
  $1.0$)

  **Unit:** $mol~m_{MP}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\gt 0$
    - **Length:** 1
  ```

`GIEX_REFQ`

: Reference solid phase concentration (optional, defaults to
  $1.0$)

  **Unit:** $mol~m_{SP}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\gt 0$
    - **Length:** 1
  ```

`GIEX_PHREFC0`

: Reference liquid phase concentration used for the pH dependent part
  of the characteristic charge (optional, defaults to `GIEX_REFC0`)

  **Unit:** $mol~m_{MP}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\gt 0$
    - **Length:** 1
  ```

`GIEX_PHREFQ`

: Reference solid phase concentration used for the pH dependent part
  of the characteristic charge (optional, defaults to `GIEX_REFQ`)

  **Unit:** $mol~m_{SP}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\gt 0$
    - **Length:** 1
  ```
