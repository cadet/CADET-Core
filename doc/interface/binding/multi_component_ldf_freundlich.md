(multi-component-ldf-freundlich-config)=

# Multi Component Linear Driving Force Freundlich

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = MULTI_COMPONENT_LDF_FREUNDLICH**

For information on model equations, refer to [](#multi-component-ldf-freundlich-model).

`IS_KINETIC`

: Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
  quasi-stationary. If a single value is given, the mode is set for all
  bound states. Otherwise, the adsorption mode is set for each bound
  state separately.

|   |   |   |
| --- | --- | --- |
| **Type:** int | **Range:** {0,1} |  |

`MCLDFFRL_KLDF`

: Rate constants in linear driving force approach

**Unit:** $s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`MCLDFFRL_KF`

: Proportionality constants

**Unit:** $m_{MP}^{3/n}~m_{SP}^{-3}~mol^{1-1/n}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`MCLDFFRL_EXP`

: Freundlich exponent

**Unit:** [-]

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\gt 0$ |  |

`MCLDFFRL_A`

: Component influences in row-major ordering

**Unit:** [-]

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ | **Length:** $\text{NCOMP}^2$ |

`MCLDFFRL_TAU`

: Small constant that ensures numerical stability

**Unit:** $mol~m_{MP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\gt 0$ |  |
