(spline-interpolation-config)=

# Spline Interpolation

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = SPLINE_INTERPOLATION**

For information on model equations, refer to [](#spline-interpolation).

`IS_KINETIC`

: Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
  quasi-stationary. If a single value is given, the mode is set for all
  bound states. Otherwise, the adsorption mode is set for each bound
  state separately.

|   |   |   |
| --- | --- | --- |
| **Type:** int | **Range:** {0,1} | **Length:** 1/NTOTALBND |

`CP_VALS_COMP_XXX`

: Pore-phase concentration support points for component `XXX` used to
  construct the spline isotherm. The values must be given in ascending
  order.

**Unit:** $mol~m_{MP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\geq 0$ |  |

`CS_VALS_COMP_XXX`

: Solid-phase equilibrium loadings corresponding to `C_VALS_COMP_XXX`
  for component `XXX`. This parameter is used if component `XXX` has
  exactly one bound state.

**Unit:** $mol~m_{SP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\geq 0$ |  |

`CS_VALS_COMP_XXX_BND_YYY`

: Solid-phase equilibrium loadings corresponding to `C_VALS_COMP_XXX`
  for bound state `YYY` of component `XXX`. This parameter is used if
  component `XXX` has more than one bound state.

**Unit:** $mol~m_{SP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\geq 0$ |  |

`SPLINE_KKIN`

: Linear-driving-force coefficients in component-major ordering.

**Unit:** $s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\geq 0$ |  |
