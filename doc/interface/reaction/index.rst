.. _FFReaction:

Reaction models
===============

Reactions can take place on different phases in the model. The interface specification depends on the phase(s) for which the reaction takes place.

A reaction can be located in:

- **Unit Operation**: The reaction group is defined under the corresponding transport phase i.e. ``unit_XXX/liquid_reaction_YYY``.
- **Particle**: The reaction group is defined under the corresponding particle phase i.e. ``particle_type_XXX/liquid_reaction_YYY`` or ``particle_type_XXX/solid_reaction_YYY``.

Furthermore, we differ between reactions which take place in either a single phase or between two phases.

Single Phase Reactions
----------------------
Single phase reactions take place either in the liquid phase or the solid phase.

The number of reaction models is specified at the unit or particle level using the following parameters:

  - ``NREAC_LIQUID`` for liquid phase reactions i.e. ``unit_XXX/NREAC_LIQUID`` or ``particle_type_XXX/NREAC_LIQUID``.,
  - ``NREAC_SOLID`` for solid phase reactions i.e. ``unit_XXX/NREAC_SOLID`` or ``particle_type_XXX/NREAC_SOLID``.

The input group for parameters of single reactions is given by ``phase_reaction_YYY`` where ``phase`` is one of ``liquid`` or ``solid`` and ``YYY`` is a 0-based index for the reaction group.

**Single-Phase Reaction Models**:

.. toctree::
    :maxdepth: 1

    mass_action_law
    michaelis_menten_kinetics


Cross Phase Reactions
---------------------
Cross phase reactions take place between a liquid phase and a solid phase.

To specify the number of reaction models, use the parameter at the column or particle level:
  - ``NREAC_CROSS_PHASE`` for cross-phase reactions i.e. ``unit_XXX/NREAC_CROSS_PHASE`` or ``particle_type_XXX/NREAC_CROSS_PHASE``.
  
The input group for cross-phase reactions is given by ``cross_phase_reaction_YYY`` where ``YYY`` is a 0-based index for the reaction group.

**Cross-Phase Reaction Models**:

.. toctree::
    :maxdepth: 1

    mass_action_law_cross_phase

.. note::
  - Cross-phase reaction models can **only** be used in the cross-phase formulation

Step by Step guide to define reaction model
-------------------------------------------
1. Select the phase(s) the reaction (liquid, solid or cross phase) will take place.

2. Choose how many reaction model types you want to define and set up the corresponding parameter in the configuration file (e.g., ``NREAC_LIQUID``, ``NREAC_SOLID``, ``NREAC_CROSS_PHASE``).
  **Note**:
    - The input group depends on whether the reaction is defined in a particle phase(s) or bulk phases.
    - It is possible to define multiple reactions of the same type.

3. Choose the appropriate reaction model for each of the reaction model type based on the kinetics of your system (e.g., Mass Action Law (:ref:`mass_action_law_config`), Michaelis-Menten (:ref:`michaelis_menten_kinetics_config`)) and write it under the corresponding input group (e.g., ``liquid_reaction_000``, ``solid_reaction_001``, ``cross_phase_reaction_002``).

4. Define the reaction parameters according to the chosen model's specifications. Refer to the specific reaction model documentation for details on required parameters.


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
