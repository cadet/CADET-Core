(extended-mobile-phase-modulator-langmuir-config)=

# Extended Mobile Phase Modulator Langmuir

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = EXTENDED_MOBILE_PHASE_MODULATOR**

For information on model equations, refer to [](#extended-mobile-phase-modulator-langmuir-model).

`IS_KINETIC`

: Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
  quasi-stationary. If a single value is given, the mode is set for all
  bound states. Otherwise, the adsorption mode is set for each bound
  state separately.

|   |   |   |
| --- | --- | --- |
| **Type:** int | **Range:** {0,1} |  |

`EMPM_COMP_MODE`

: Determines the mode of each component ($0$ denotes the modifier
  component, $1$ is linear binding, $2$ is modified Langmuir
  binding). At most one modifier component is allowed, that is, a
  modifier is not required.

  Note that this field has the same name for the externally dependent
  variant of the model.

|   |   |   |
| --- | --- | --- |
| **Type:** int | **Range:** $\{0,1,2\}$ |  |

`EMPM_KA`

: Adsorption rate constants

**Unit:** $m_{MP}^3~mol^{-1}~s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`EMPM_KD`

: Desorption rate constants

**Unit:** $m_{MP}^{3\beta}~mol^{-\beta}~s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`EMPM_QMAX`

: Maximum adsorption capacities

**Unit:** $mol~m_{SP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`EMPM_BETA`

: Parameters describing the ion-exchange characteristics (IEX)

**Unit:** $mol~m_{SP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`EMPM_GAMMA`

: Parameters describing the hydrophobicity (HIC)

**Unit:** $m_{MP}^{3} mol^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\mathbb{R}$ |  |
