(return)=
# Return data

## Group /input/return

`WRITE_SOLUTION_TIMES`

: Write times at which a solution was produced (optional, defaults to 1)

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_IDAS_META`

: Write IDAS time integrator meta data (optional, defaults to 0), e.g. time step sizes

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SOLUTION_LAST`

: Write full solution state vector at last time point, including the state derivative vector (optional, defaults to 0)

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SENS_LAST`

: Write full sensitivity state vectors at last time point, including the state derivative vector (optional, defaults to 0)

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`SPLIT_COMPONENTS_DATA`

: Determines whether a joint dataset (matrix or tensor) for all components is created or if each component is put in a separate dataset ($\texttt{XXX_COMP_000}$, $\texttt{XXX_COMP_001}$, etc.) (optional, defaults to 1)

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`SPLIT_PORTS_DATA`

: Determines whether a joint dataset (matrix or tensor) for all inlet/outlet ports is created or if each port is put in a separate dataset ($\texttt{XXX_PORT_000}$, $\texttt{XXX_PORT_001}$, etc.) (optional, defaults to 1)

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`SINGLE_AS_MULTI_PORT`

: Determines whether single port unit operations are treated as multi port unit operations in the output naming scheme (i.e., $\texttt{_PORT_XYZ_}$ is added to the name) (optional, defaults to 0)

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

## Group /input/return/unit_XXX

All fields of this group are optional and default to $0$, i.e. nothing is written unless it is requested.

`WRITE_COORDINATES`

: Write coordinates of discretization nodes

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SOLUTION_INLET`

: Write solutions at unit operation inlet $c^l_i(t,0)$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SOLUTION_OUTLET`

: Write solutions at unit operation outlet (chromatograms) $c^l_i(t,L)$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SOLUTION_BULK`

: Write solutions of the bulk volume $c^l_i$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SOLUTION_PARTICLE`

: Write solutions of the particle mobile phase $c^p_{j,i}$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SOLUTION_SOLID`

: Write solutions of the solid phase $c^s_{j,i,m_{j,i}}$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SOLUTION_FLUX`

: Write solutions of the bead fluxes $j_{f,i}$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SOLUTION_VOLUME`

: Write solutions of the volume V

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SOLDOT_INLET`

: Write solution time derivatives at unit operation inlet $\partial c^l_i(t,0) / \partial t$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SOLDOT_OUTLET`

: Write solution time derivatives at unit operation outlet (chromatograms) $\partial c^l_i(t,L) / \partial t$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SOLDOT_BULK`

: Write solution time derivatives of the bulk volume $\partial c^l_i / \partial t$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SOLDOT_PARTICLE`

: Write solution time derivatives of the particle mobile phase $\partial c^p_{j,i} / \partial t$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SOLDOT_SOLID`

: Write solution time derivatives of the solid phase $\partial c^s_{j,i,m_{j,i}} / \partial t$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SOLDOT_FLUX`

: Write solution time derivatives of the bead fluxes $\partial j_{f,i} / \partial t$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SOLDOT_VOLUME`

: Write solution time derivatives of the volume $\partial V / \partial t$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SENS_INLET`

: Write sensitivities at unit operation inlet $\partial c^l_i(t,0) / \partial p$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SENS_OUTLET`

: Write sensitivities at unit operation outlet (chromatograms) $\partial c^l_i(t,L) / \partial p$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SENS_BULK`

: Write sensitivities of the bulk volume $\partial c^l_i / \partial p$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SENS_PARTICLE`

: Write sensitivities of the particle mobile phase $\partial c^p_{j,i} / \partial p$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SENS_SOLID`

: Write sensitivities of the solid phase $\partial c^s_{j,i,m_{j,i}} / \partial p$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SENS_FLUX`

: Write sensitivities of the bead fluxes $\partial j_{f,i} / \partial p$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SENS_VOLUME`

: Write sensitivities of the volume $\partial V / \partial p$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SENSDOT_INLET`

: Write sensitivity time derivatives at unit operation inlet $\partial^2 c^l_i(t,0) / (\partial p, \partial t)$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SENSDOT_OUTLET`

: Write sensitivity time derivatives at unit operation outlet (chromatograms) $\partial^2 c^l_i(t,L) / (\partial p, \partial t)$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SENSDOT_BULK`

: Write sensitivity time derivatives of the bulk volume $\partial^2 c^l_i / (\partial p, \partial t)$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SENSDOT_PARTICLE`

: Write sensitivity time derivatives of the particle mobile phase $\partial^2 c^p_{j,i} / (\partial p, \partial t)$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SENSDOT_SOLID`

: Write sensitivity time derivatives of the solid phase $\partial^2 c^s_{j,i,m_{j,i}} / (\partial p, \partial t)$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SENSDOT_FLUX`

: Write sensitivity time derivatives of the bead fluxes $\partial^2 j_{f,i} / (\partial p, \partial t)$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SENSDOT_VOLUME`

: Write sensitivity time derivatives of the volume $\partial^2 V / (\partial p, \partial t)$

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SOLUTION_LAST_UNIT`

: Write solution state vector of this unit at last time point, including the state derivative vector (optional, defaults to 0)

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```

`WRITE_SENS_LAST_UNIT`

: Write sensitivity state vector of this unit at last time point, including the state derivative vector (optional, defaults to 0)

  ```{list-table}

  * - **Type:** int
    - **Range:** $\{0,1\}$
  ```
