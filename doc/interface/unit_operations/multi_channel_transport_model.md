(multi-channel-transport-model-config)=

# Multichannel Transport model (MCT model)

## Group /input/model/unit_XXX - UNIT_TYPE = MULTI_CHANNEL_TRANSPORT

For information on model equations, refer to [](#multi-channel-transport-model-model).

`UNIT_TYPE`

> Specifies the type of unit operation model
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** string | **Range:** `MULTI_CHANNEL_TRANSPORT` | **Length:** 1 |

`NCOMP`

> Number of chemical components in the chromatographic medium
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 1$ | **Length:** 1 |

`NREAC_LIQUID`

> Number of liquid phase reaction models (optional, only if liquid reactions are present).
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 0$ | **Length:** 1 |

## Group /input/model/unit_XXX/liquid_reaction_YYY

Each reaction is specified in another subgroup `liquid_reaction_YYY`, see [](#FFReaction).

`INIT_C`

> Initial concentrations for each component in all channels of the bulk mobile phase (length `NCOMP`), or for each component in each channel (length $\texttt{NCOMP} \cdot \texttt{NCHANNEL}$, ordering channel-major)
>
> **Unit:** $\mathrm{mol}\,\mathrm{m}^{-3}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** :math:`\texttt{NCOMP} / \texttt{NCOMP} \cdot \texttt{NCHANNEL |

`INIT\_STATE`

> Full state vector for initialization (optional, `INIT\_C` will be ignored; if length is $2\texttt{NDOF}$, then the second half is used for time derivatives).
> The ordering of the state vector is defined in [](#UnitOperationStateOrdering), similar to 2D models but with channels instead of radial column direction.
>
> **Unit:** $various$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** $\texttt{NDOF} / 2\texttt{NDOF}$ |

`COL_DISPERSION`

> Axial dispersion coefficient. In case of a spatially inhomogeneous setting, the `SENS_PARTYPE` field is used for indexing the channel when specifying parameter sensitivities.
>
> **Unit:** $\mathrm{m}^{2}\,\mathrm{s}^{-1}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** see `COL_DISPERSION_MULTIPLEX` |

`COL_DISPERSION_MULTIPLEX`

> > Multiplexing mode of `COL_DISPERSION`. Determines whether `COL_DISPERSION` is treated as component-, channel-, and/or section-independent. This field is optional. When left out, multiplexing behavior is inferred from the length of `COL_DISPERSION`. Valid modes are:
>
> 0. Component-independent, channel-independent, section-independent; length of `COL_DISPERSION` is 1
> 1. Component-independent, channel-dependent, section-independent; length of `COL_DISPERSION` is `NCHANNEL`
> 2. Component-dependent, channel-independent, section-independent; length of `COL_DISPERSION` is `NCOMP`
> 3. Component-dependent, channel-dependent, section-independent; length of `COL_DISPERSION` is $\texttt{NCOMP} \cdot \texttt{NCHANNEL}$; ordering is channel-major
> 4. Component-independent, channel-independent, section-dependent; length of `COL_DISPERSION` is `NSEC`
> 5. Component-independent, channel-dependent, section-dependent; length of `COL_DISPERSION` is $\texttt{NCHANNEL} \cdot \texttt{NSEC}$; ordering is section-major
> 6. Component-dependent, channel-independent, section-independent; length of `COL_DISPERSION` is $\texttt{NCOMP} \cdot \texttt{NSEC}$; ordering is section-major
> 7. Component-dependent, channel-dependent, section-dependent; length of `COL_DISPERSION` is $\texttt{NCOMP} \cdot \texttt{NCHANNEL} \cdot \texttt{NSEC}$; ordering is section-channel-major
>
> > |   |   |   |
> > | --- | --- | --- |
> > | **Type:** int | **Range:** $\{0, \dots, 7 \}$ | **Length:** 1 |

`COL_LENGTH`

> Column length
>
> **Unit:** $\mathrm{m}$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $> 0$ | **Length:** 1 |

`VELOCITY`

> Velocity
>
> **Unit:** $\mathrm{m}\,\mathrm{s}^{-1}$
>
> Indicates flow direction in each channel (forward if value is positive, backward if value is negative), see Section [](#MUOPGRMflow2D)). In case of a spatially inhomogeneous setting, the `SENS_PARTYPE` field is used for indexing the channel when specifying parameter sensitivities.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\mathbb{R}$ | **Length:** see `VELOCITY_MULTIPLEX` |

`VELOCITY_MULTIPLEX`

> > Multiplexing mode of `VELOCITY`. Determines whether `VELOCITY` is treated as channel- and/or section-independent. This field is optional. When left out, multiplexing behavior is inferred from the length of `VELOCITY`. Valid modes are:
>
> 0. Channel-independent, section-independent; length of `VELOCITY` is 1
> 1. Channel-dependent, section-independent; length of `VELOCITY` is `NCHANNEL`
> 2. Section-dependent; length of `VELOCITY` is `NSEC`
> 3. Channel-dependent, section-dependent; length of `VELOCITY` is $\texttt{NCHANNEL} \cdot \texttt{NSEC}$; ordering is section-major
>
> > |   |   |   |
> > | --- | --- | --- |
> > | **Type:** int | **Range:** $\{0, \dots, 3 \}$ | **Length:** 1 |

`EXCHANGE_TYPE`

> Type of exchange model to be used (optional)
>
> **Unit:** N/A
>
> Specifies the exchange model to be applied for inter-channel transport. Possible values are:
>
> - `LINEAR_EX`: Linear Exchange model (default)
> - `LANGMUIR_EX`: Langmuir Exchange model
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** string | **Range:** {LANGMUIR_EX, LINEAR_EX} | **Length:** 1 |

`EXCHANGE_MATRIX`

> > Exchange matrix
> >
> > **Unit:** $\mathrm{s}^{-1}$
> >
> > Ordered list containing all exchange rates $e^k_{ij}$ for component $k$ from compartment $i$ to $j$ based on the exchange matrix $E^k$. The vector ordering is source channel - destination channel - component (i.e. i-j-k) major.
> >
> > $$
> > E^k=\begin{bmatrix}
> > 0 & e^k_{12} & \dots & e^k_{1N} \\
> > e^k_{21} & \ddots & & \vdots\\
> > \vdots & & \ddots & e^k_{(N-1)N}\\
> > e^k_{N1} & \dots & e^k_{N(N-1)} & 0
> > \end{bmatrix}
> > $$
> >
> > For addressing the exchange rates as a parameter senstivity, the mapping is as follows:
>
> - `SENS_BOUNDPHASE` *Channel from*
> - `SENS_PARTYPE` *Channel to*
>
> > |   |   |   |
> > | --- | --- | --- |
> > | **Type:** double | **Range:** $[0,1]$ | **Length:** :math:`\texttt{NCHANNEL}*\texttt{NCHANNEL}*\texttt{NCOMP} |

`SATURATION_MATRIX`

: Maximum exchange saturation matrix (optional)
  Only when using the Langmuir Exchange model.

  **Unit:** $mol~m_{channel}^{-3}$

  Matrix containing maximum saturation levels $q_{max,i}^{k}$ for each component $i$ in each channel $k$.
  This parameter determines the saturation limit for the Langmuir-type exchange kinetics.

  - If $q_{max,i}^{k} = 0$ the corresponding exchange flux will reduce to the linear exchange model.
  - If $q_{max,i}^{k} < 0$ the model will describe an Anti-Langmuir behavior.

  |   |   |   |
  | --- | --- | --- |
  | **Type:** double | **Range:** $\mathbb{R}$ |  |

`NCHANNEL`

> Number of channels $ij$
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 1$ | **Length:** 1 |

`CHANNEL_CROSS_SECTION_AREAS`

> Cross section areas
>
> **Unit:** $\mathrm{m}^{2}$
>
> Defines the cross section area of each channel. The `SENS_PARTYPE` field is used for indexing the channel when specifying parameter sensitivities.
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $> 0$ | **Length:** `NCHANNEL` |

## Group /input/model/unit_XXX/discretization - UNIT_TYPE = MULTI_CHANNEL_TRANSPORT

`USE_ANALYTIC_JACOBIAN`

> Determines whether analytically computed Jacobian matrix (faster) is used (value is 1) instead of Jacobians generated by algorithmic differentiation (slower, value is 0)
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\{0, 1\}$ | **Length:** 1 |

`NCOL`

> Number of axial column discretization cells
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** int | **Range:** $\geq 1$ | **Length:** 1 |

`LINEAR_SOLVER_BULK`

> > Linear solver used for the sparse column bulk block. This field is optional, the best available method is selected (i.e., sparse direct solver if possible). Valid values are:
>
> - `DENSE` Converts the sparse matrix into a banded matrix and uses regular LAPACK. Slow and memory intensive, but always available.
> - `UMFPACK` Uses the UMFPACK sparse direct solver (LU decomposition) from SuiteSparse. Fast, but has to be enabled when compiling and requires UMFPACK library.
> - `SUPERLU` Uses the SuperLU sparse direct solver (LU decomposition). Fast, but has to be enabled when compiling and requires SuperLU library.
>
> > |   |   |   |
> > | --- | --- | --- |
> > | **Type:** string | **Range:** $\{\texttt{DENSE},\texttt{UMFPACK},\texttt{SUPERLU}\}$ | **Length:** 1 |

`RECONSTRUCTION`

> Type of reconstruction method for FV fluxes
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** string | **Range:** `WENO, KOREN` | **Length:** 1 |

For further discretization parameters (including `GRID_FACES` for non-equidistant grids), see also [](#flux-reconstruction-methods), and [](#non-consistency-solver-parameters).
