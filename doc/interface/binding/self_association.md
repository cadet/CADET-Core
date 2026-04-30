(self-association-config)=

# Self Association

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = SELF_ASSOCIATION**

For information on model equations, refer to [](#self-association-model).

`IS_KINETIC`

: Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
  quasi-stationary. If a single value is given, the mode is set for all
  bound states. Otherwise, the adsorption mode is set for each bound
  state separately.

|   |   |   |
| --- | --- | --- |
| **Type:** int | **Range:** {0,1} |  |

`SAI_KA1`

: Adsorption rate constants

**Unit:** $m_{MP}^3~m_{SP}^{-3}~s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`SAI_KA2`

: Adsorption rate constants

**Unit:** $m_{MP}^6~m_{SP}^{-3}~mol^{-1}~s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`SAI_KD`

: Desorption rate constants

**Unit:** $s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`SAI_NU`

: Characteristic charges $\nu$ of the protein

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`SAI_SIGMA`

: Steric factors $\sigma$ of the protein

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`SAI_LAMBDA`

: Stationary phase capacity (monovalent salt counterions); The total
  number of binding sites available on the resin surface

**Unit:** $mol m_{SP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`SAI_REFC0`

: Reference liquid phase concentration (optional, defaults to
  $1.0$)

**Unit:** $mol m_{MP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\gt 0$ |  |

`SAI_REFQ`

: Reference solid phase concentration (optional, defaults to
  $1.0$)

**Unit:** $mol m_{SP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\gt 0$ |  |
