(ffmodelsystem)=
# System of unit operations

## Group /input/model

`INIT_STATE_Y`

: Initial full state vector (optional, unit operation specific initial data is ignored)

  ```{list-table}

  * - **Type:** double
    - **Length:** $\texttt{NDOF}$
  ```

`INIT_STATE_YDOT`

: Initial full time derivative state vector (optional, unit operation specific initial data is ignored)

  ```{list-table}

  * - **Type:** double
    - **Length:** $\texttt{NDOF}$
  ```

`INIT_STATE_SENSY_XXX`

: Initial full state vector of the $\texttt{XXX}$ th sensitivity system (optional, can currently not be specified on unit operation level)

  ```{list-table}

  * - **Type:** double
    - **Length:** $\texttt{NDOF}$
  ```

`INIT_STATE_SENSYDOT_XXX`

: Initial full time derivative state vector of the $\texttt{XXX}$ th sensitivity system (optional, can currently not be specified on unit operation level)

  ```{list-table}

  * - **Type:** double
    - **Length:** $\texttt{NDOF}$
  ```

(ffmodelsystemconnections)=
## Group /input/model/connections

`NSWITCHES`

: Number of valve switches

  ```{list-table}

  * - **Type:** int
    - **Range:** $\geq 1$
    - **Length:** 1
  ```

`CONNECTIONS_INCLUDE_PORTS`

: Determines whether the $\texttt{CONNECTIONS}$ table includes ports ($1$) or not ($0$). Optional, defaults to 0 unless a unit operation model with multiple ports is present.

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{ 0,1 \}$
    - **Length:** 1
  ```

`CONNECTIONS_INCLUDE_DYNAMIC_FLOW`

: Determines whether the $\texttt{CONNECTIONS}$ table includes linear, quadratic, and cubic flow rate coefficients (1) or not (0). Optional, defaults to 0.

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{ 0,1 \}$
    - **Length:** 1
  ```

(ffmodelconnectionswitch)=
## Group /input/model/connections/switch_XXX

`SECTION`

: Index of the section that activates this connection set

  ```{list-table}

  * - **Type:** int
    - **Range:** $\geq 0$
    - **Length:** 1
  ```

`CONNECTIONS`

: Matrix with list of connections in row-major storage. Columns are *UnitOpID from*, *UnitOpID to*, *Port from*, *Port to*, *Component from*, *Component to*, *volumetric flow rate*, *linear flow rate coefficient*, *quadratic flow rate coefficient*, *cubic flow rate coefficient*.
  If both port indices are $-1$, all ports are connected.
  If both component indices are $-1$, all components are connected.

  The flow rate is a cubic function of time,

  ```{math}
  Q = Q_0 + Q_1(t - t_s) + Q_2(t-t_s)^2 + Q_3(t-t_s)^3,
  ```

  where $t_s$ is the beginning of the section that activates the switch (i.e., $\texttt{SECTION_TIMES}$ at index $\texttt{SECTION}$).

  The port indices are left out if $\texttt{CONNECTIONS_INCLUDE_PORTS}$ is set to $0$ and no unit operation with multiple ports is present in the system. If a unit operation with multiple ports is present, $\texttt{CONNECTIONS_INCLUDE_PORTS}$ is ignored and port indices are mandatory.

  The last three flow rate coefficients are left out if $\texttt{CONNECTIONS_INCLUDE_DYNAMIC_FLOW}$ is set to $0$.
  Contrary to the constant coefficient, which has the parameter name $\texttt{CONNECTION}$, the other coefficients are named $\texttt{CONNECTION_LIN}$, $\texttt{CONNECTION_QUAD}$, and $\texttt{CONNECTION_CUB}$, respectively.

  For addressing the flow rates as a parameter sensitivity, the mapping is as follows:

  - $\texttt{SENS_UNIT}$ Unused, always set to $-1$
  - $\texttt{SENS_BOUNDPHASE}$ *UnitOpID from*
  - $\texttt{SENS_REACTION}$ *UnitOpID to*
  - $\texttt{SENS_COMP}$ *Port from*
  - $\texttt{SENS_PARTYPE}$ *Port to*
  - $\texttt{SENS_SECTION}$ $\texttt{SECTION}$ that activates the valve switch

  ```{list-table}

  * - **Type:** double
    - **Range:** $\geq -1$
    - **Length:** $\{5,7,8,10\} \cdot \texttt{NCONNECTIONS}$
  ```

(ffmodelexternalsourcelininterp)=
## Group /input/model/external/source_XXX - EXTFUN_TYPE = LINEAR_INTERP_DATA

`VELOCITY`

: Velocity of the external profile in positive column axial direction.
  The velocity is normalized to a column with length 1, hence the unit $\mathrm{s}^{-1}$.
  A negative value lets the profile travel towards the column inlet instead.
  The value must not be $0$, since the profile is evaluated at the transformed time
  $(1 - z) / \texttt{VELOCITY} + t$.

  **Unit:** $\mathrm{s}^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R} \setminus \{0\}$
    - **Length:** 1
  ```

`DATA`

: Function values $T$ at the data points

  **Unit:** $[\mathrm{Ext}]$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** Arbitrary
  ```

`TIME`

: Time of the data points

  **Unit:** $\mathrm{s}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\geq 0.0$
    - **Length:** Same as $\texttt{DATA}$
  ```

(ffmodelexternalsourcepiececubicpoly)=
## Group /input/model/external/source_XXX - EXTFUN_TYPE = PIECEWISE_CUBIC_POLY

`VELOCITY`

: Velocity of the external profile in positive column axial direction.
  The velocity is normalized to a column with length 1, hence the unit $\mathrm{s}^{-1}$.
  A negative value lets the profile travel towards the column inlet instead.
  The value must not be $0$, since the profile is evaluated at the transformed time
  $(1 - z) / \texttt{VELOCITY} + t$.

  **Unit:** $\mathrm{s}^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R} \setminus \{0\}$
    - **Length:** 1
  ```

`CONST_COEFF`

: Constant coefficients of piecewise cubic polynomial

  **Unit:** $[\mathrm{Ext}]$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** Arbitrary
  ```

`LIN_COEFF`

: Linear coefficients of piecewise cubic polynomial

  **Unit:** $[\mathrm{Ext}]\,\mathrm{s}^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** Same as $\texttt{CONST_COEFF}$
  ```

`QUAD_COEFF`

: Quadratic coefficients of piecewise cubic polynomial

  **Unit:** $[\mathrm{Ext}]\,\mathrm{s}^{-2}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** Same as $\texttt{CONST_COEFF}$
  ```

`CUBE_COEFF`

: Cubic coefficients of piecewise cubic polynomial

  **Unit:** $[\mathrm{Ext}]\,\mathrm{s}^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** Same as $\texttt{CONST_COEFF}$
  ```

`SECTION_TIMES`

: Simulation times at which a new piece begins (breaks of the piecewise polynomial)

  **Unit:** $\mathrm{s}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\geq 0.0$
    - **Length:** $\texttt{CONST_COEFF}+1$
  ```

(ffmodelsolver)=
## Group /input/model/solver

`GS_TYPE`

: Type of Gram-Schmidt orthogonalization, see IDAS guide Section 4.5.7.3, p. 41f. A value of $0$ enables classical Gram-Schmidt, a value of 1 uses modified Gram-Schmidt.

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0, 1\}$
    - **Length:** 1
  ```

`MAX_KRYLOV`

: Defines the size of the Krylov subspace in the iterative linear GMRES solver that solves the
  Schur complement of the network. That Schur complement is formed on the inlet degrees of freedom
  of the network, of which there are
  $\texttt{NCOUPLINGDOF} = \sum_i \texttt{NPORT}_i \cdot \texttt{NCOMP}_i$, summed over all
  unit operations that have an inlet (0: $\texttt{MAX_KRYLOV} = \texttt{NCOUPLINGDOF}$).

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0, \dots, \texttt{NCOUPLINGDOF}\}$
    - **Length:** 1
  ```

`MAX_RESTARTS`

: Maximum number of restarts in the GMRES algorithm. If lack of memory is not an issue, better use a larger Krylov space than restarts.

  ```{list-table}

  * - **Type:** int
    - **Range:** $\geq 0$
    - **Length:** 1
  ```

`SCHUR_SAFETY`

: Schur safety factor; Influences the tradeoff between linear iterations and nonlinear error control; see IDAS guide Section 2.1 and 5.

  ```{list-table}

  * - **Type:** double
    - **Range:** $\geq 0$
    - **Length:** 1
  ```

`LINEAR_SOLUTION_MODE`

: Determines whether the system of models is solved in parallel (1) or sequentially (2). A sequential solution is only possible for systems without cyclic connections. The setting can be chosen automatically (0) based on a heuristic (less than 25 unit operations and acyclic network selects sequential mode). Optional, defaults to automatic (0).

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{ 0,1,2 \}$
    - **Length:** 1
  ```
