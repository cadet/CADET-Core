.. _FFReaction:

Reaction models
===============

.. toctree::
    :maxdepth: 2

    mass_action_law
    michaelis_menten_kinetics

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


.. _multiple-particle-types_reactions:

Multiple particle types
-----------------------

The group that contains the parameters of a reaction model in unit operation with index ``XXX`` reads ``/input/model/unit_XXX/reaction_particle``.
This is valid for models with a single particle type.
If a model has multiple particle types, it may have a different reaction model in each type.
The parameters are then placed in the group ``/input/model/unit_XXX/reaction_particle_YYY`` instead, where ``YYY`` denotes the index of the particle type.

Note that, in any case, ``/input/model/unit_XXX/reaction_particle_000`` contains the parameters of the first (and possibly sole) particle type.
This group also takes precedence over a possibly existing ``/input/model/unit_XXX/adsorption_particle`` group.
