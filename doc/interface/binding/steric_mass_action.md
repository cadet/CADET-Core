(steric-mass-action-config)=

# Steric Mass Action

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = STERIC_MASS_ACTION**

For information on model equations, refer to [](#steric-mass-action-model).

`IS_KINETIC`

: Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
  quasi-stationary. If a single value is given, the mode is set for all
  bound states. Otherwise, the adsorption mode is set for each bound
  state separately.

|   |   |   |
| --- | --- | --- |
| **Type:** int | **Range:** {0,1} |  |

`SMA_KA`

: Adsorption rate constants

**Unit:** $m_{MP}^3~m_{SP}^{-3}~s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`SMA_KD`

: Desorption rate constants

**Unit:** $s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`SMA_NU`

: Characteristic charges of the protein; The number of sites
  $\nu$ that the protein interacts with on the resin surface

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`SMA_SIGMA`

: Steric factors of the protein; The number of sites $\sigma$ on
  the surface that are shielded by the protein and prevented from
  exchange with the salt counterions in solution

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`SMA_LAMBDA`

: Stationary phase capacity (monovalent salt counterions); The total
  number of binding sites available on the resin surface

**Unit:** $mol~m_{SP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`SMA_REFC0`

: Reference liquid phase concentration (optional, defaults to
  $1.0$)

**Unit:** $mol~m_{MP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\gt 0$ |  |

`SMA_REFQ`

: Reference solid phase concentration (optional, defaults to
  $1.0$)

**Unit:** $mol~m_{SP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\gt 0$ |  |
