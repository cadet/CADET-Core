(bi-steric-mass-action-config)=

# Bi Steric Mass Action

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = BI_STERIC_MASS_ACTION**

For information on model equations, refer to [](#bi-steric-mass-action-model).

`IS_KINETIC`

: Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
  quasi-stationary. If a single value is given, the mode is set for all
  bound states. Otherwise, the adsorption mode is set for each bound
  state separately.

|   |   |   |
| --- | --- | --- |
| **Type:** int | **Range:** {0,1} | **Length:** 1/NTOTALBND |

`BISMA_KA`

: Adsorption rate constants in state-major ordering

**Unit:** $m_{MP}^{3}~m_{SP}^{-3}~s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`BISMA_KD`

: Desorption rate constants in state-major ordering

**Unit:** $s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`BISMA_NU`

: Characteristic charges $\nu_{i,j}$ of the $i$th protein
  with respect to the $j$th binding site type in state-major
  ordering

**Unit:** $s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`BISMA_SIGMA`

: Steric factors $\sigma_{i,j}$ of the $i$th protein with
  respect to the $j$th binding site type in state-major
  ordering

**Unit:** $s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`BISMA_LAMBDA`

: Stationary phase capacity (monovalent salt counterions) of the
  different binding site types $\lambda_j$

**Unit:** $mol~m_{SP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`BISMA_REFC0`

: Reference liquid phase concentration for each binding site type or
  one value for all types (optional, defaults to $1.0$)

**Unit:** $mol~m_{MP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\gt 0$ |  |

`BISMA_REFQ`

: Reference solid phase concentration for each binding site type or one
  value for all types (optional, defaults to $1.0$)

**Unit:** $mol~m_{SP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\gt 0$ |  |
