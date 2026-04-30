(inlet-config)=

# Inlet

## Group /input/model/unit_XXX - UNIT-TYPE = INLET

For information on model equations, refer to [](#inlet-model).

`UNIT_TYPE`

> Specifies the type of unit operation model
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** string | **Range:** `INLET` | **Length:** 1 |

`NCOMP`

> Number of chemical components in the chromatographic medium
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 1$ | **Length:** 1 |

`INLET_TYPE`

> Specifies the type of inlet profile
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** string | **Range:** `PIECEWISE_CUBIC_POLY` | **Length:** 1 |

## Group /input/model/unit_XXX/sec_XXX

`CONST_COEFF`

> Constant coefficients for inlet concentrations
>
> **Unit:** $\mathrm{mol}\,\mathrm{m}_{\mathrm{IV}}^{-3}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** `NCOMP` |

`LIN_COEFF`

> Linear coefficients for inlet concentrations
>
> **Unit:** $\mathrm{mol}\,\mathrm{m}_{\mathrm{IV}}^{-3}\,\mathrm{s}^{-1}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** `NCOMP` |

`QUAD_COEFF`

> Quadratic coefficients for inlet concentrations
>
> **Unit:** $\mathrm{mol}\,\mathrm{m}_{\mathrm{IV}}^{-3}\,\mathrm{s}^{-2}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** `NCOMP` |

`CUBE_COEFF`

> Cubic coefficients for inlet concentrations
>
> **Unit:** $\mathrm{mol}\,\mathrm{m}_{\mathrm{IV}}^{-3}\,\mathrm{s}^{-3}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** `NCOMP` |
