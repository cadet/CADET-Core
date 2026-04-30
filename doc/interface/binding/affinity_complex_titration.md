(affinity-complex-titration-config)=

# Affinity Complex Titration

## General remarks

- The **first component** is the mobile phase modulator.
- This first component must be **non-binding**.
- The ACT implementation currently supports **at most one bound state per component**.
- The first component can be specified either as a negative logarithmic concentration
  ($\mathrm{pH}$, $\mathrm{pNa}$, ...) or as a raw ion concentration.
- To convert an existing ACT setup from `ACT_USE_ION_CONC = False` to `ACT_USE_ION_CONC = True` while keeping exactly the same model response:
  : 1. Replace the first component value $\mathrm{pIon}$ by $c_{\mathrm{ion}} = 10^{-\mathrm{pIon}}$.
    2. Replace `ACT_PKAA` by `ACT_CMID_A = 10^{-\mathrm{ACT\_PKAA}}`.
    3. Replace `ACT_PKAG` by `ACT_CMID_G = 10^{-\mathrm{ACT\_PKAG}}`.
    4. Keep `ACT_KA`, `ACT_KD`, `ACT_QMAX`, `ACT_ETAA`, and `ACT_ETAG` unchanged.

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = AFFINITY_COMPLEX_TITRATION**

For information on model equations, refer to [](#affinity-complex-titration).

`IS_KINETIC`

: Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 = quasi-stationary.
  If a single value is given, the mode is set for all bound states. Otherwise, the adsorption
  mode is set for each bound state separately.

|   |   |   |
| --- | --- | --- |
| **Type:** int | **Range:** {0,1} |  |

`ACT_USE_ION_CONC`

: Selects how the mobile phase modulator concentration is interpreted.

  - `False`: the first component concentration is on a negative logarithmic axis, for example
    $\mathrm{pH}$ or $\mathrm{pNa}$.
    Use `ACT_PKAA` and `ACT_PKAG`.
  - `True`: the first component is the raw ion concentration $c_{\mathrm{ion}}$.
    Use `ACT_CMID_A` and `ACT_CMID_G`.

  Default is `False`.

  These two cases are equivalent when

  $$
  \mathrm{p}K_{a,A,i} = -\log_{10}(c_{\mathrm{mid},A,i}), \qquad
  \mathrm{p}K_{a,G,i} = -\log_{10}(c_{\mathrm{mid},G,i}).
  $$

|   |   |   |
| --- | --- | --- |
| **Type:** bool | **Range:** {False, True} |  |

`ACT_KA`

: Adsorption rate constants.

**Unit:** $m_{MP}^3\,mol^{-1}\,s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`ACT_KD`

: Desorption rate constants.

**Unit:** $s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`ACT_QMAX`

: Maximum binding capacities before modulation by the ACT capacity gate.

**Unit:** $mol\,m_{SP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $> 0$ |  |

`ACT_ETAA`

: Hill-type coefficients controlling how strongly the apparent binding capacity changes with the modulator concentration.

**Unit:** $1$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\mathbb{R}$ |  |

`ACT_ETAG`

: Hill-type coefficients controlling how strongly the apparent affinity changes with the modulator concentration.

**Unit:** $1$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\mathbb{R}$ |  |

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ACT_USE_ION_CONC = False**

Use these parameters only when the first component is given as $\mathrm{pIon}$
(for example $\mathrm{pH}$).

`ACT_PKAA`

: Midpoint of the binding capacity transition on the negative logarithmic ion axis.

**Unit:** $1$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\mathbb{R}$ |  |

`ACT_PKAG`

: Midpoint of the affinity (equilibrium constant) transition on the negative logarithmic ion axis.

**Unit:** $1$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\mathbb{R}$ |  |

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ACT_USE_ION_CONC = True**

Use these parameters only when the first component is given as a raw ion concentration.
CADET internally converts them to the same negative logarithmic axis used by the `ACT_PKAA` / `ACT_PKAG` form.

`ACT_CMID_A`

: Midpoint ion concentration for the binding capacity transition.
  Must be non-negative and should use the same concentration unit as the first liquid-phase component.

**Recommended unit:** $mol\,m_{MP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`ACT_CMID_G`

: Midpoint ion concentration for the affinity (equilibrium constant) transition.
  Must be non-negative and should use the same concentration unit as the first liquid-phase component.

**Recommended unit:** $mol\,m_{MP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |
