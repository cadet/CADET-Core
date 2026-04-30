(hic-water-on-hydrophobic-surfaces-config)=

# HIC Water on Hydrophobic Surfaces

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = HIC_WATER_ON_HYDROPHOBIC_SURFACES**

For information on model equations, refer to [](#hic-water-on-hydrophobic-surfaces-model).

`IS_KINETIC`

: Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
  quasi-stationary. If a single value is given, the mode is set for all
  bound states. Otherwise, the adsorption mode is set for each bound
  state separately.

|   |   |   |
| --- | --- | --- |
| **Type:** int | **Range:** {0,1} | **Length:** 1/NTOTALBND |

`HICWHS_KA`

: Adsorption rate constant

**Unit:** $m_{MP}^{3}~m_{SP}^{-3}~s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`HICWHS_KD`

: Desorption rate constant

**Unit:** $s^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`HICWHS_NU`

: Number of ligands per ligand-protein interaction

**Unit:** [-]

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`HICWHS_QMAX`

: Maximum binding capacity

**Unit:** $mol~m_{SP}^{-3}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`HICWHS_BETA0`

: Parameters describing the number of highly ordered water molecules
  that stabilize the hydrophobic surfaces at infinitely diluted
  salt concentration

**Unit:** [-]

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |

`HICWHS_BETA1`

: Parameters describing the change in the number of highly ordered
  water molecules that stabilize the hydrophobic surfaces with
  respect to changes in the salt concentration

**Unit:** $m_{MP}^{3}~mol^{-1}$

|   |   |   |
| --- | --- | --- |
| **Type:** double | **Range:** $\ge 0$ |  |
