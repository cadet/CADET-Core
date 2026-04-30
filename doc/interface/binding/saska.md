(saska-config)=

# Saska

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = SASKA**

For information on model equations, refer to [](#saska-model).

`IS_KINETIC`

: Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
  quasi-stationary. If a single value is given, the mode is set for all
  bound states. Otherwise, the adsorption mode is set for each bound
  state separately.

|   |   |   |
| --- | --- | --- |
| **Type:** int | **Range:** {0,1} |  |

`SASKA_H`

: Henry coefficient

**Unit:** $m_{MP}^3~m_{SP}^{-3}~s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\mathbb {R}$ |  |

`SASKA_K`

: Quadratic factors

**Unit:** $m_{MP}^6~m_{SP}^{-3}~s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\mathbb {R}$ |  |
