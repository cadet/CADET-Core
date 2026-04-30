(pbm-config)=

# Crystallization / Precipitation models

Crystallization / Precipitation in CADET can be modeled by [](#primary-particle-formation), [](#aggregation), [](#fragmentation), as well as their combinations.
The configuration is described for all three modules in the following.
For more information on the model equations, please refer to [](#FFCrystallization).

Note that any of these models can be used in any unit operation that supports reactions.

:::{attention}
By setting the `NCOMP` field, you specify the number of FV cells for the internal coordinate:
If primary particle formation is considered, we have two components that account for the solute $c$ and solubility $c_\text{eq}$.
Note that the first component must be solute $c$ and the last component must be the solubility $c_\text{eq}$.
Additionally, we discretize the particle size domain (internal coordinate) using customary FV methods, giving us a finite set of particle sizes under consideration $\{x_1, \dots, x_{N_x}\}$.
Every particle size considered is treated as an individual component of the unit operation and the field `NCOMP` of that unit operation in which the crystallization happens, must be specified accordingly as $N_x$ or $N_x + 2$ if the model includes primary particle formation.
:::

Example code for configuring the crystallization models is available in [CADET-Verification](https://github.com/cadet/CADET-Verification/) .

## Group /input/model/unit_XXX

`NCOMP`

> Number of components, which is defined by the number of considered particle sizes (i.e. discretization with FV cells) plus, in the case that primary particle formation is employed, 2 (for solute and solubility).
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 3$ | **Length:** 1 |

`NREAC_LIQUID`

> Number of liquid phase reaction models (optional, only if liquid reactions are present).
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 0$ | **Length:** 1 |

## Group /input/model/unit_XXX/liquid_reaction_000 - TYPE = CRYSTALLIZATION - UNIT_TYPE = CSTR

*The following parameters need to be specified under Group /input/model/unit_XXX/liquid_reaction_000/*

`CRY_MODE`

> > Crystallization mode, which determines the exact model equation to be employed.
>
> 1. Pure primary particle formation, as described in [](#primary-particle-formation).
> 2. Pure aggregation, as described in [](#aggregation).
> 3. Combined primary particle formation and aggregation.
> 4. Pure fragmentation, as described in [](#fragmentation).
> 5. Combined primary particle formation and fragmentation.
> 6. Combined aggregation and fragmentation.
> 7. Combined primary particle formation, aggregation, and fragmentation.
>
> > |   |   |   |
> > | --- | --- | --- |
> > | **Type:** double | **Range:** $\geq 1$ | **Length:** $\mathrm{N_x} + 1$ |

`CRY_BINS`

> Coordinates of the cell faces, e.g. equidistant or logarithmic discretization of the internal coordinate $x \in [x_c, x_\text{end}]$, including the end points.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 1$ | **Length:** $\mathrm{N_x} + 1$ |

### Population and Mass Balance input parameters

`CRY_NUCLEI_MASS_DENSITY`

> Nulcei mass density
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** 1 |

`CRY_VOL_SHAPE_FACTOR`

> Volumetric shape factor of the particles
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** 1 |

`CRY_PRIMARY_NUCLEATION_RATE`

> Primary nucleation rate constant $k_p$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** 1 |

`CRY_SECONDARY_NUCLEATION_RATE`

> Secondary nucleation rate $k_b$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** 1 |

`CRY_GROWTH_RATE_CONSTANT`

> Growth rate constant $k_g$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** 1 |

`CRY_GROWTH_CONSTANT`

> Growth constant $\gamma$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** 1 |

`CRY_A`

> Defines constant $a$ used to determine the growth rate
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** 1 |

`CRY_G`

> Defines constant $g$ used to determine the growth rate
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** 1 |

`CRY_P`

> Defines constant $p$ used to determine the growth rate
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** 1 |

`CRY_GROWTH_DISPERSION_RATE`

> Growth dispersion rate $D_g$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** 1 |

`CRY_U`

> Defines constant $u$ used to determine the primary nucleation
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** 1 |

`CRY_B`

> Defines constant $b$ used to determine the secondary nucleation
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** 1 |

`CRY_K`

> Defines constant $k$ used to determine the secondary nucleation, usually set to $\geq 1$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** 1 |

`CRY_GROWTH_SCHEME_ORDER`

> Defines the growth flux FV reconstruction scheme. It can only be
>
> - $1$: upwind scheme
> - $2$: HR Koren scheme
> - $3$: WENO23 scheme
> - $4$: WENO35 scheme.
>
> We recommend using the HR Koren scheme, which showed to be the most performant in our benchmarks.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $[1, \dots, 4]$ | **Length:** 1 |

### Aggregation input parameters

`CRY_AGGREGATION_INDEX`

> Defines the aggregation kernel. It can only be
>
> - $0$: constant kernel
> - $1$: Brownian kernel
> - $2$: Smoluchowski kernel
> - $3$: Golovin kernel
> - $4$: differential force kernel
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $[0, \dots, 4]$ | **Length:** 1 |

`CRY_AGGREGATION_RATE_CONSTANT`

> Aggregation rate constant $\beta_0$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $> 0$ | **Length:** 1 |

### Fragmentation input parameters

`CRY_FRAGMENTATION_RATE_CONSTANT`

> Fragmentation rate constant $S_0$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $> 0$ | **Length:** 1 |

`CRY_FRAGMENTATION_KERNEL_GAMMA`

> Fragmentation kernel coefficient $\gamma$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $> 1.0$ | **Length:** 1 |

`CRY_FRAGMENTATION_SELECTION_FUNCTION_ALPHA`

> Fragmentation selection function coefficient $\alpha$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $> 0$ | **Length:** 1 |
