.. _FFCrystallization:

Crystallization / Precipitation Models
======================================

CADET features various crystallization / precipitation mechanisms:

 - :ref:`pbm_model`
 - :ref:`aggregation_model`
 - :ref:`fragmentation_model`

 Additionally, any combination of these models is featured.
 Further, any of these crystallization \ precipitation models can be combined with any of the transport models listed among the :ref:`unit_operation_models`.
 To this end, most common use-cases consider a :ref:`cstr_model` or a dispersive plug-flow reactor (which can be modeled with the :ref:`lumped_rate_model_without_pores_model`).

 For detailed information on the crystallization models implemented in CADET, including nucleation, growth, aggregation and fragmentation, please refer to :cite:`Zhang2024` and :cite:`Zhang2025`.
