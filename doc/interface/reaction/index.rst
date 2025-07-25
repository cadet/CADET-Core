.. _FFReaction:

Reaction models
===============

This section describes the configuration of reaction models for different phases in CADET unit operations.

The corresponding group that contains the reaction models for each phase reads:

- For bulk phase: ``/input/model/unit_XXX/reaction_bulk``
- For solid phase: ``/input/model/unit_XXX/reaction_solid``
- For particle phase: ``/input/model/unit_XXX/reaction_particle_YYY``
- For cross phase: ``/input/model/unit_XXX/reaction_cross_phase_YYY``

Where ``XXX`` denotes the index of the unit operation and ``YYY`` denotes the index of the particle type.

Configuration Steps
-------------------

Before defining the reaction model parameters, it is necessary to specify how many reaction models are used in each phase.

1. **Set the number of reactions**: This is done by setting ``/NREACT`` to the number of reaction models used in the phase.

   For example:

   - ``/input/model/unit_000/reaction_bulk/NREACT = 2`` indicates that 2 reaction models are used in the bulk phase
   - ``/input/model/unit_000/reaction_solid/NREACT = 1`` indicates that 1 reaction model is used in the solid phase

2. **Define reaction model parameters**: After setting ``NREACT``, the parameters for each individual reaction model are specified.

The group that contains the parameters of each reaction model reads:
``/input/model/unit_XXX/reaction_phase/reaction_model_ZZZ``. 

Where ``XXX`` denotes the index of the unit operation, ``phase`` is either ``bulk``, ``solid``, ``particle_YYY``,
or ``cross_phase_YYY``, and ``ZZZ`` denotes the index of the reaction model (starting from ``000``).

**Example**: If ``NREACT = 2`` for the bulk phase, then the parameters would be found in:

- ``/input/model/unit_000/reaction_bulk/reaction_model_000`` (first reaction model)
- ``/input/model/unit_000/reaction_bulk/reaction_model_001`` (second reaction model)

The reaction model parameters are specified in the following sections:

**Single-Phase Reaction Models** (for bulk, solid, and particle phases):

.. toctree::
    :maxdepth: 2

    mass_action_law
    michaelis_menten_kinetics

**Cross-Phase Reaction Models** (only for cross-phase reactions):

.. toctree::
    :maxdepth: 2

    mass_action_law_cross_phase

.. note::
   Cross-phase reaction models can **only** be used in the cross-phase configuration (``reaction_cross_phase_YYY``).
   They cannot be used in bulk, solid, or particle phase configurations.
   For reactions within a single phase, use the standard reaction models listed above.


Externally dependent reaction models
------------------------------------

Some reaction models have a variant that can use external sources as specified :ref:`/input/model/external/<FFModelExternalSourceLinInterp>` (also see Section :ref:`dependence-on-external-function_react`).
For the sake of brevity, only the standard variant of those reaction models is specified below.
In order to obtain the format for the externally dependent variant, first replace the reaction model name ``XXX`` by ``EXT_XXX``.
Each parameter :math:`p` (except for stoichiometric and exponent matrices) depends on a (possibly distinct) external source in a polynomial way:

.. math::

    \begin{aligned}
        p(T) &= p_{\texttt{TTT}} T^3 + p_{\texttt{TT}} T^2 + p_{\texttt{T}} T + p.
    \end{aligned}

Thus, a parameter ``XXX_YYY`` of the standard reaction model variant is replaced by the four parameters ``EXT_XXX_YYY``, ``EXT_XXX_YYY_T``, ``EXT_XXX_YYY_TT``, and ``EXT_XXX_YYY_TTT``.
Since each parameter can depend on a different external source, the dataset ``EXTFUN`` (not listed in the standard variants below) should contain a vector of 0-based integer indices of the external source of each parameter.
The ordering of the parameters in ``EXTFUN`` is given by the ordering in the standard variant.
However, if only one index is passed in ``EXTFUN``, this external source is used for all parameters.

Note that parameter sensitivities with respect to column radius, column length, particle core radius, and particle radius may be wrong when using externally dependent reaction models.
This is caused by not taking into account the derivative of the external profile with respect to column position.
