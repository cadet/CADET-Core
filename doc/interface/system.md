(ffmodelsystem)=

# System of unit operations

## Group /input/model

`NUNITS`

> Number of unit operations in the system
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 1$ | **Length:** 1 |

`INIT_STATE_Y`

> Initial full state vector (optional, unit operation specific initial data is ignored)
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 1$ | **Length:** 1 |

`INIT_STATE_YDOT`

> Initial full time derivative state vector (optional, unit operation specific initial data is ignored)
>
> |   |   |
> | --- | --- |
> | **Type:** double |  |

`INIT_STATE_SENSY_XXX`

> Initial full time derivative state vector of the `XXX` th sensitivity system (optional, can currently not be specified on unit operation level)
>
> |   |   |
> | --- | --- |
> | **Type:** double |  |

`INIT_STATE_SENSYDOT_XXX`

> Initial full state vector of the `XXX` th sensitivity system (optional, can currently not be specified on unit operation level)
>
> |   |   |
> | --- | --- |
> | **Type:** double |  |

(ffmodelsystemconnections)=

## Group /input/model/connections

`NSWITCHES`

> Number of valve switches
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 1$ | **Length:** 1 |

`CONNECTIONS_INCLUDE_PORTS`

> Determines whether the `CONNECTIONS` table includes ports ($1$) or not ($0$). Optional, defaults to 0 unless a unit operation model with multiple ports is present.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\{ 0,1 \}$ | **Length:** 1 |

`CONNECTIONS_INCLUDE_DYNAMIC_FLOW`

> Determines whether the `CONNECTIONS` table includes linear, quadratic, and cubic flow rate coefficients (1) or not (0). Optional, defaults to 0.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\{ 0,1 \}$ | **Length:** 1 |

(ffmodelconnectionswitch)=

## Group /input/model/connections/switch_XXX

`SECTION`

> Index of the section that activates this connection set
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 0$ | **Length:** 1 |

`CONNECTIONS`

> > Matrix with list of connections in row-major storage. Columns are *UnitOpID from*, *UnitOpID to*, *Port from*, *Port to*, *Component from*, *Component to*, *volumetric flow rate*, *linear flow rate coefficient*, *quadratic flow rate coefficient*, *cubic flow rate coefficient*.
> > If both port indices are $-1$, all ports are connected.
> > If both component indices are $-1$, all components are connected.
> >
> > The flow rate is a cubic function of time,
> >
> > $$
> > Q = Q_0 + Q_1(t - t_s) + Q_2(t-t_s)^2 + Q_3(t-t_s)^3,
> > $$
> >
> > where $t_s$ is the beginning of the section that activates the switch (i.e., `SECTION_TIMES` at index `SECTION`).
> >
> > The port indices are left out if `CONNECTIONS_INCLUDE_PORTS` is set to $0$ and no unit operation with multiple ports is present in the system. If a unit operation with multiple ports is present, `CONNECTIONS_INCLUDE_PORTS` is ignored and port indices are mandatory.
> >
> > The last three flow rate coefficients are left out if `CONNECTIONS_INCLUDE_DYNAMIC_FLOW` is set to $0$.
> > Contrary to the constant coefficient, which has the parameter name `CONNECTION`, the other coefficients are named `CONNECTION_LIN`, `CONNECTION_QUAD`, and `CONNECTION_CUB`, respectively.
> >
> > For addressing the flow rates as a parameter senstivity, the mapping is as follows:
>
> - `SENS_UNIT` Unused, always set to $-1$
> - `SENS_BOUNDPHASE` *UnitOpID from*
> - `SENS_REACTION` *UnitOpID to*
> - `SENS_COMP` *Port from*
> - `SENS_PARTYPE` *Port to*
> - `SENS_SECTION` `SECTION` that activates the valve switch
>
> > |   |   |   |
> > | --- | --- | --- |
> > | **Type:** double | **Range:** $\geq -1$ | **Length:** $\{5,7,8,10\} \cdot \texttt{NCONNECTIONS}$ |

(ffmodelexternalsourcelininterp)=

## Group /input/model/external/source_XXX - EXTFUN_TYPE = LINEAR_INTERP_DATA

`VELOCITY`

> Velocity of the external profile in positive column axial direction.
> The velocity is normalized to a column with length 1, hence the unit $\mathrm{s}^{-1}$.
>
> **Unit:** $\mathrm{s}^{-1}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** 1 |

`DATA`

> Function values $T$ at the data points
>
> **Unit:** $[\mathrm{Ext}]$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** Arbitrary |

`TIME`

: Time of the data points

  **Unit:** $\mathrm{s}$

  |   |   |   |
  | --- | --- | --- |
  | **Type:** double | **Range:** $\geq 0.0$ | **Length:** Same as `DATA` |

(ffmodelexternalsourcepiececubicpoly)=

## Group /input/model/external/source_XXX - EXTFUN_TYPE = PIECEWISE_CUBIC_POLY

`VELOCITY`

> Velocity of the external profile in positive column axial direction.
> The velocity is normalized to a column with length 1, hence the unit $\mathrm{s}^{-1}$.
>
> **Unit:** $\mathrm{s}^{-1}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** 1 |

`CONST_COEFF`

> Constant coefficients of piecewise cubic polynomial
>
> **Unit:** $[\mathrm{Ext}]$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** Arbitrary |

`LIN_COEFF`

> Linear coefficients of piecewise cubic polynomial
>
> **Unit:** $[\mathrm{Ext}]\,\mathrm{s}^{-1}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** Same as `CONST_COEFF` |

`QUAD_COEFF`

> Quadratic coefficients of piecewise cubic polynomial
>
> **Unit:** $[\mathrm{Ext}]\,\mathrm{s}^{-2}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** Same as `CONST_COEFF` |

`CUBE_COEFF`

> Cubic coefficients of piecewise cubic polynomial
>
> **Unit:** $[\mathrm{Ext}]\,\mathrm{s}^{-3}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** Same as `CONST_COEFF` |

`SECTION_TIMES`

> Simulation times at which a new piece begins (breaks of the piecewise polynomial)
>
> **Unit:** $\mathrm{s}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0.0$ | **Length:** $\texttt{CONST\_COEFF}+1$ |

(ffmodelsolver)=

## Group /input/model/solver

`GS_TYPE`

> Type of Gram-Schmidt orthogonalization, see IDAS guide Section~4.5.7.3, p.~41f. A value of $0$ enables classical Gram-Schmidt, a value of 1 uses modified Gram-Schmidt.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\{0, 1\}$ | **Length:** 1 |

`MAX_KRYLOV`

> Defines the size of the Krylov subspace in the iterative linear GMRES solver (0: $\texttt{MAX\_KRYLOV} = \texttt{NDOF}$)
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\{0, \dots, \texttt{NDOF}\}$ | **Length:** 1 |

`MAX_RESTARTS`

> Maximum number of restarts in the GMRES algorithm. If lack of memory is not an issue, better use a larger Krylov space than restarts.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 0$ | **Length:** 1 |

`SCHUR_SAFETY`

> Schur safety factor; Influences the tradeoff between linear iterations and nonlinear error control; see IDAS guide Section~2.1 and 5.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** 1 |

`LINEAR_SOLUTION_MODE`

> Determines whether the system of models is solved in parallel (1) or sequentially (2). A sequential solution is only possible for systems without cyclic connections. The setting can be chosen automatically (0) based on a heuristic (less than 25 unit operations and acyclic network selects sequential mode). Optional, defaults to automatic (0).
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\{ 0,1,2 \}$ | **Length:** 1 |
