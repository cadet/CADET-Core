(sensitivity)=

# Parameter Sensitivities

(ffsensitivity)=

## Group /input/sensitivity

`NSENS`

> Number of sensitivities to be computed
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 0$ | **Length:** 1 |

`SENS_METHOD`

> Method used for computation of sensitivities (algorithmic differentiation)
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** string | **Range:** `ad1` | **Length:** 1 |

(ffsensitivityparam)=

## Group /input/sensitivity/param_XXX

`SENS_UNIT`

> Unit operation index
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 0$ | **Length:** $\geq 1$ |

`SENS_NAME`

> Name of the parameter (Note that `PAR_RADIUS` and `PAR_CORERADIUS` sensitivities are only available for Finite Volume discretization)
>
> |   |   |
> | --- | --- |
> | **Type:** string |  |

`SENS_COMP`

> Component index ($-1$ if parameter is independent of components, see the multiplexing of the corresponding parameter, e.g. `FILM_DIFFUSION_MULTIPLEX`)
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq -1$ |  |

`SENS_PARTYPE`

> Particle type index ($-1$ if parameter is independent of particle types, see the multiplexing of the corresponding parameter, e.g. `PORE_DIFFUSION_MULTIPLEX`)
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq -1$ |  |

`SENS_REACTION`

> Reaction index ($-1$ if parameter is independent of reactions)
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq -1$ |  |

`SENS_BOUNDPHASE`

> Bound phase index ($-1$ if parameter is independent of bound phases)
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq -1$ | **Length:** $\geq 1$ |

`SENS_SECTION`

> Section index ($-1$ if parameter is independent of sections)
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq -1$ | **Length:** $\geq 1$ |

`SENS_ABSTOL`

> Absolute tolerance used in the computation of the sensitivities (optional). Rule of thumb: $\texttt{ABSTOL} / \texttt{PARAM\_VALUE}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0.0$ | **Length:** $\geq 1$ |

`SENS_FACTOR`

> Linear factor of the combined sensitivity (optional, taken as $1.0$ if left out)
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** $\geq 1$ |
