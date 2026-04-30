(multi-component-anti-langmuir-config)=

# Multi Component Anti-Langmuir

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = MULTI_COMPONENT_ANTILANGMUIR**

For information on model equations, refer to [](#multi-component-anti-langmuir-model).

`IS_KINETIC`

: Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
  quasi-stationary. If a single value is given, the mode is set for all
  bound states. Otherwise, the adsorption mode is set for each bound
  state separately.

|   |   |   |
| --- | --- | --- |
| **Type:** int | **Range:** {0,1} |  |

`MCAL_KA`

: Adsorption rate constants

**Unit:** $m_{MP}^3~mol^{-1}~s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`MCAL_KD`

: Desorption rate constants

**Unit:** $s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`MCAL_QMAX`

: Maximum adsorption capacities

**Unit:** $mol~m_{SP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\gt 0$ |  |

`MCAL_ANTILANGMUIR`

: Anti-Langmuir coefficients (optional)

**Unit:** $mol~m_{SP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** {-1,1} |  |
