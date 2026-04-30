(ffcrystallization)=

# Crystallization / Precipitation Models

CADET-Core features various crystallization / precipitation mechanisms, modeled using population balance models (PBM).
These models are defined by a particle-number continuity equation which describes the evolution of the number density $n$ of the particles over time $t$ and with respect to size, the so-called internal coordinate $x$, and external coordinate $z$.
The external coordinate can be a characteristic dimension of the reactor itself, including its axial length.

The PBM in CADET-Core are implemented in a modular manner and can be used in any unit operation that supports reactions.
Typical applications consider crystallization in a CSTR or, to model continuous processes, in a Dispersive Plug-Flow Reactor (DPFR), which can be described by the LRM without solid phase.

For detailed information on the PBM implemented in CADET-Core, please refer to {cite}`Zhang2024` and {cite}`Zhang2025`.

Additionally, any combination of these models is also supported.
Further, any of these crystallization / precipitation models can be combined with any of the transport models listed among the [unit operation models](#unit-operation-models).
To this end, most common use-cases consider a [CSTR](#cstr-model) or a dispersive plug-flow reactor (which can be modeled with the [LRM without pores](#lumped-rate-model-without-pores-model)).

```{tableofcontents}
```
