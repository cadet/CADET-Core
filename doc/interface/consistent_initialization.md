(non-consistency-solver-parameters)=
# Nonlinear solver for consistent initialization

## Group /input/model/unit_XXX/discretization/consistency_solver - Nonlinear consistency solver parameters

The following specifications of the nonlinear solver for consistent initialization are available.
The default is a `COMPOSITE` solver that applies multiple solvers subsequently, where the default is `ATRN_ERR` followed by `LEVMAR`.

`SOLVER_NAME`

: Name of the solver.
  Available solvers are

  - `LEVMAR`: Levenberg-Marquardt method (ported from the MATLAB implementation)

  - `ATRN_RES`: Adaptive trust-region Newton solver, implementation of the NLEQ-RES algorithm described in {cite}`Deuflhard2011` (p. 131)

  - `ATRN_ERR`: Robust adaptive trust-region Newton solver, implementation of the NLEQ-ERR algorithm described in {cite}`Deuflhard2011` (pp. 148)

  - `PTC`: Pseudo-transient continuation method, see {cite}`Deuflhard2011` (section 6.4)

  - `COMPOSITE`: Applies multiple solvers subsequently

  ```{list-table}

  * - **Type:** string
    - **Length:** $1$
  ```

`SUBSOLVERS`

: Vector with names of solvers for the `COMPOSITE` solver (only required for composite solver). See `SOLVER_NAME` for available solvers.

  ```{list-table}

  * - **Type:** string
    - **Length:** $\gt 1$
  ```

`INIT_DAMPING`

: Initial damping factor (default is $0.01$), only required for `LEVMAR`, `ATRN_RES` and `ATRN_ERR`. See `SOLVER_NAME` for available solvers.

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** $1$
  ```

`MIN_DAMPING`

: Minimal damping factor (default is $0.0001$), only required for `ATRN_RES` and `ATRN_ERR`. See `SOLVER_NAME` for available solvers.

  ```{list-table}

  * - **Type:** double
    - **Range:** $\ge 0$
    - **Length:** $1$
  ```

`MAX_ITERATIONS`

: Maximum number of iterations (optional, defaults to $50$ for `LEVMAR`, `ATRN_RES` and `ATRN_ERR`, and to $100$ for `PTC`).

  ```{list-table}

  * - **Type:** int
    - **Range:** $\ge 1$
    - **Length:** $1$
  ```

`INIT_PSEUDOTIMESTEP`

: Initial pseudo time step (optional, defaults to $20.0$), only used by `PTC`.

  ```{list-table}

  * - **Type:** double
    - **Range:** $> 0$
    - **Length:** $1$
  ```

`MAX_NONMONOTONE_ITERATIONS`

: Maximum number of consecutive iterations in which the residual is allowed to grow (optional, defaults to $5$), only used by `PTC`.

  ```{list-table}

  * - **Type:** int
    - **Range:** $\ge 0$
    - **Length:** $1$
  ```

`NUMERIC_VARIANT`

: Selects the numeric variant of the method (optional, defaults to $0$), only used by `PTC`.

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
    - **Length:** $1$
  ```
