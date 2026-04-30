(parameter-dependencies)=

# Parameter Dependencies

Some parameters depend on other parameters (parameter-parameter dependency) or the solution variables (parameter-state dependency).
Parameter dependencies are defined in the unit operation scope.

## Parameter-Parameter Dependencies

### Group /input/model/unit_XXX

`COL_DISPERSION_DEP`

> Parameter dependence of column dispersion on the interstitial velocity. Available for the LRM, LRMP and GRM units (with FV discretization only at the moment)
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** string | **Range:** `POWER_LAW` | **Length:** 1 |

### Group /input/model/unit_XXX/particle_type_YYY

`FILM_DIFFUSION_DEP`

> Parameter dependence of film diffusion on the interstitial velocity. Available for the LRMP unit (with FV discretization only at the moment)
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** string | **Range:** `POWER_LAW` | **Length:** 1 |

#### **Correlations**

Different types of parameter correlations are can be applied.
The following correlations can be used for all parameter-parameter dependencies, but we specify the required input fields only for `COL_DISPERSION_DEP`, for the sake of conciseness.

**Power Law**

$$
\begin{aligned}
    p_{dep} &= p_{dep} \cdot b \ |p_{on}^x|
\end{aligned}
$$

Here, $p_{dep}$ is the dependent parameter and $p_{on}$ is the parameter it depends on.

`COL_DISPERSION_DEP_BASE`

> Base $b$ of the power law parameter dependence. Optional, defaults to $1.0$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |

`COL_DISPERSION_DEP_EXPONENT`

> Exponent $x$ of the power law parameter dependence
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |

`COL_DISPERSION_DEP_ABS`

> Specifies whether or not the absolute value should be computed. Optional, defaults to $1$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\{0, 1\}$ | **Length:** 1 |

## Parameter-State Dependencies

Currently, the dependence of surface diffusion on the particle liquid salt component is the only available parameter-state dependence.

### Group /input/model/unit_XXX/particle_type_YYY

`SURFACE_DIFFUSION_DEP`

> Parameter dependence of `SURFACE_DIFFUSION` on the particle liquid salt component (i.e. component with index 0). Valid dependencies are:
>
> - `NONE` Original parameter is used unmodified.
> - `LIQUID_SALT_EXPONENTIAL` Original parameter is modified by exponential law of liquid phase salt concentration.
> - `LIQUID_SALT_POWER` Original parameter is modified by power law of liquid phase salt concentration.
> - `LIQUID_SALT_COLLOIDAL_AFFINITY` Original parameter is modified by colloidal binding affinity based on liquid phase salt concentration.
>
> Optional: If left out, no parameter dependence is assumed and the original surface diffusion coefficients are used unmodified.
>
> |   |   |
> | --- | --- |
> | **Type:** string | **Length:** $1 / \texttt{NPARTYPE}$ |

`SURFACE_DIFFUSION_EXPFACTOR`

> Factor `p1` in exponential law particle surface diffusion relation
> $D_{s, i, m} = \tilde{D}_{s, i, m} p_{1, i, m} exp \left(p_{2, i, m} c_{0}^{p} \right)$, where $\tilde{D}_{s, i, m}$ is the original surface diffusion coefficient.
> Only required if `SURFACE_DIFFUSION_DEP` is `LIQUID_SALT_EXPONENTIAL`.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** `NBOUND` |
>
> `SURFACE_DIFFUSION_EXPFACTOR`
> $D_{s, i, m} = \tilde{D}_{s, i, m} \left[  p_{4, i, m} \left( k_{i, m} \left( c_{0}^{p} \right) \right)^{p_{5, i, m}} p_{6, i, m} exp \left( p_{7, i, m} k_{i, m} \left( c_{0}^{p} \right) \right) \right]$
> where $\tilde{D}_{s, i, m}$ is the original surface diffusion coefficient and
> $k_{i, m} \left( c_{0}^{p} \right) = p_{1, i, m}\left( c_{0}^{p} \right)^{p_{2, i, m}} + p_{3, i, m}$.
> Only required if `SURFACE_DIFFUSION_DEP` is `LIQUID_SALT_COLLOIDAL_AFFINITY`.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** `NBOUND` |

`SURFACE_DIFFUSION_EXPARGMULT`

> Factor `p2` in exponential law particle surface diffusion relation
> $D_{s, i, m} = \tilde{D}_{s, i, m} p_{1, i, m} exp \left(p_{2, i, m} c_{0}^{p} \right)$
> where $\tilde{D}_{s, i, m}$ is the original surface diffusion coefficient. Only required if `SURFACE_DIFFUSION_DEP` is `LIQUID_SALT_EXPONENTIAL`.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** `NBOUND` |

`SURFACE_DIFFUSION_POWFACTOR`

> Factor `p1` in power law particle surface diffusion relation
> $D_{s, i, m} = \tilde{D}_{s, i, m} p_{1, i, m} \left( c_{0}^{p} \right)^{p_{2, i, m}}$
> where $\tilde{D}_{s, i, m}$ is the original surface diffusion coefficient. Only required if `SURFACE_DIFFUSION_DEP` is `LIQUID_SALT_POWER`.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** `NBOUND` |

`SURFACE_DIFFUSION_POWEXP`

> Fjactor `p2` in power law particle surface diffusion relation
> $D_{s, i, m} = \tilde{D}_{s, i, m} p_{1, i, m} \left( c_{0}^{p} \right)^{p_{2, i, m}}$
> where $\tilde{D}_{s, i, m}$ is the original surface diffusion coefficient. Only required if `SURFACE_DIFFUSION_DEP` is `LIQUID_SALT_POWER`.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** `NBOUND` |

`SURFACE_DIFFUSION_LOGKEQFACTOR`

> Factor `p1` in colloidal affinity law particle surface diffusion relation
> $D_{s, i, m} = \tilde{D}_{s, i, m} \left[  p_{4, i, m} \left( k_{i, m} \left( c_{0}^{p} \right) \right)^{p_{5, i, m}} p_{6, i, m} exp \left( p_{7, i, m} k_{i, m} \left( c_{0}^{p} \right) \right) \right]$
> where $\tilde{D}_{s, i, m}$ is the original surface diffusion coefficient and
> $k_{i, m} \left( c_{0}^{p} \right) = p_{1, i, m}\left( c_{0}^{p} \right)^{p_{2, i, m}} + p_{3, i, m}$.
> Only required if `SURFACE_DIFFUSION_DEP` is `LIQUID_SALT_COLLOIDAL_AFFINITY`.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** `NBOUND` |

`SURFACE_DIFFUSION_LOGKEQEXP`

> Factor `p2` in colloidal affinity law particle surface diffusion relation
> $D_{s, i, m} = \tilde{D}_{s, i, m} \left[  p_{4, i, m} \left( k_{i, m} \left( c_{0}^{p} \right) \right)^{p_{5, i, m}} p_{6, i, m} exp \left( p_{7, i, m} k_{i, m} \left( c_{0}^{p} \right) \right) \right]$
> where $\tilde{D}_{s, i, m}$ is the original surface diffusion coefficient and
> $k_{i, m} \left( c_{0}^{p} \right) = p_{1, i, m}\left( c_{0}^{p} \right)^{p_{2, i, m}} + p_{3, i, m}$.
> Only required if `SURFACE_DIFFUSION_DEP` is `LIQUID_SALT_COLLOIDAL_AFFINITY`.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** `NBOUND` |

`SURFACE_DIFFUSION_LOGKEQCONST`

> Factor `p3` in colloidal affinity law particle surface diffusion relation
> $D_{s, i, m} = \tilde{D}_{s, i, m} \left[  p_{4, i, m} \left( k_{i, m} \left( c_{0}^{p} \right) \right)^{p_{5, i, m}} p_{6, i, m} exp \left( p_{7, i, m} k_{i, m} \left( c_{0}^{p} \right) \right) \right]$
> where $\tilde{D}_{s, i, m}$ is the original surface diffusion coefficient and
> $k_{i, m} \left( c_{0}^{p} \right) = p_{1, i, m}\left( c_{0}^{p} \right)^{p_{2, i, m}} + p_{3, i, m}$.
> Only required if `SURFACE_DIFFUSION_DEP` is `LIQUID_SALT_COLLOIDAL_AFFINITY`.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** `NBOUND` |
