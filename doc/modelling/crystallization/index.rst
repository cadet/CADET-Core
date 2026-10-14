.. _FFCrystallization:

Crystallization / Precipitation Models
======================================

CADET-Core features various crystallization / precipitation mechanisms, modeled using population balance models (PBM). These models are defined by a particle-number continuity equation which describes the evolution of the number density :math:`n` of the particles over time :math:`t` and with respect to size, the so-called internal coordinate :math:`x`, and external coordinate :math:`z`.
The external coordinate can be a characteristic dimension of the reactor itself, including its axial length.

The PBM in CADET-Core are implemented in a modular manner and can be used in any unit operation that supports reactions.
Typical applications consider crystallization in a CSTR or, to model continuous processes, in a Dispersive Plug-Flow Reactor (DPFR), which can be described by the LRM without solid phase.

For detailed information on the PBM implemented in CADET-Core, please refer to :cite:`Zhang2024` and :cite:`Zhang2025`.
A concise overview of the models is provided in the following sections:

 .. toctree::
    :maxdepth: 1

    primary_particle_formation
    aggregation
    fragmentation

 Additionally, any combination of these models is also supported.
 Further, any of these crystallization \ precipitation models can be combined with any of the transport models listed among the :ref:`unit_operation_models`.
 To this end, most common use-cases consider a :ref:`cstr_model` or a dispersive plug-flow reactor (which can be modeled with the :ref:`lumped_rate_model_without_pores_model`).

 For detailed information on the crystallization models implemented in CADET, including PBM, aggregation and fragmentation, please refer to :cite:`Zhang2024` and :cite:`Zhang2025`.
