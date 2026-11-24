(inlet-config)=
# Inlet

## Group /input/model/unit_XXX - UNIT_TYPE = INLET

For information on model equations, refer to {ref}`inlet-model`.

`UNIT_TYPE`

: Specifies the type of unit operation model

  ```{list-table}

  * - **Type:** string
    - **Range:** $\texttt{INLET}$
    - **Length:** 1
  ```

`NCOMP`

: Number of chemical components in the chromatographic medium

  ```{list-table}

  * - **Type:** int
    - **Range:** $\geq 1$
    - **Length:** 1
  ```

`INLET_TYPE`

: Specifies the type of inlet profile

  ```{list-table}

  * - **Type:** string
    - **Range:** $\texttt{PIECEWISE_CUBIC_POLY}$
    - **Length:** 1
  ```

## Group /input/model/unit_XXX/sec_XXX

`CONST_COEFF`

: Constant coefficients for inlet concentrations

  **Unit:** $\mathrm{mol}\,\mathrm{m}_{\mathrm{IV}}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** $\texttt{NCOMP}$
  ```

`LIN_COEFF`

: Linear coefficients for inlet concentrations

  **Unit:** $\mathrm{mol}\,\mathrm{m}_{\mathrm{IV}}^{-3}\,\mathrm{s}^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** $\texttt{NCOMP}$
  ```

`QUAD_COEFF`

: Quadratic coefficients for inlet concentrations

  **Unit:** $\mathrm{mol}\,\mathrm{m}_{\mathrm{IV}}^{-3}\,\mathrm{s}^{-2}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** $\texttt{NCOMP}$
  ```

`CUBE_COEFF`

: Cubic coefficients for inlet concentrations

  **Unit:** $\mathrm{mol}\,\mathrm{m}_{\mathrm{IV}}^{-3}\,\mathrm{s}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** $\texttt{NCOMP}$
  ```
