(mass-action-law-config)=

# Mass Action Law

**Group /input/model/unit_XXX(/particle_type_YYY)/phase_reaction_ZZZ/ - TYPE = MASS_ACTION_LAW**

For information on model equations, refer to [](#mass-action-law-model).

## Notes

- `phase_reaction_ZZZ` refers to one of the phase-specific raction groups listed in [](#FFReaction), e.g., `liquid_reaction_ZZZ`, `solid_reaction_ZZZ`.
- Some dimensions below depend on the hosting phase of this model instance:
  : - Bulk phase or particle liquid phase: `NVAR = NCOMP`
    - Particle solid phase: `NVAR = NTOTALBOUND` (total number of bound states across all components)

`MAL_KFWD`

> Forward rate constants for reactions in this phase (available for external functions)
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** `NREACT` |

`MAL_KBWD`

> Backward rate constants for reactions in this phase (available for external functions)
>
> |   |   |   |
> | --- | --- | --- |
> | **Type:** double | **Range:** $\geq 0$ | **Length:** `NREACT` |

`MAL_STOICHIOMETRY`

> Stoichiometric matrix as $\texttt{NVAR} \times \texttt{NREACT}$ matrix in row-major storage (phase-dependent `NVAR`, see Notes)
>
> |   |   |
> | --- | --- |
> | **Type:** double |  |

`MAL_EXPONENTS_FWD`

> Forward exponent matrix as $\texttt{NVAR} \times \texttt{NREACT}$ matrix in row-major storage (optional, calculated from `MAL_STOICHIOMETRY` by default)
>
> |   |   |
> | --- | --- |
> | **Type:** double |  |

`MAL_EXPONENTS_BWD`

> Backward exponent matrix as $\texttt{NVAR} \times \texttt{NREACT}$ matrix in row-major storage (optional, calculated from `MAL_STOICHIOMETRY` by default)
>
> |   |   |
> | --- | --- |
> | **Type:** double |  |

## CADET Python Interface Example

This example shows the configuration of one mass-action-law reaction system.
The system has two components A and B, where A reactions to B.
In addition to that the model includes:
\- A forward reaction rate constant of 1.0
\- A backward reaction rate constant of 1.0

```

```

\# Reaction in a bulk liquid phase:

```
input.model.unit_000.NREAC_LIQUID = 1
input.model.unit_000.reaction_liquid_000.type = MASS_ACTION_LAW
input.model.unit_000.reaction_liquid_000.MAL_KFWD = [1.0]
input.model.unit_000.reaction_liquid_000.MAL_KBWD = [1.0]
input.model.unit_000.reaction_liquid_000.MAL_STOICHIOMETRY = [-1.0, 1.0] # A -> B
```

\# Reaction in a particle liquid phase:

```
input.model.unit_000.particle_type_000.NREAC_PORE = 1
input.model.unit_000.particle_type_000.reaction_liquid_000.type = MASS_ACTION_LAW
input.model.unit_000.particle_type_000.reaction_liquid_000.MAL_KFWD = [1.0]
input.model.unit_000.particle_type_000.reaction_liquid_000.MAL_KBWD = [1.0]
input.model.unit_000.particle_type_000.reaction_liquid_000.MAL_STOICHIOMETRY = [-1.0, 1.0] # A -> B
```

## See also

- Cross-phase variant: [](#mass-action-law-cross-phase-config) (supports liquid/solid modifiers and bulk coupling)
