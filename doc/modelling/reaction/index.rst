.. _reaction_models:

Reaction models
===============


Reaction models describe the (net) fluxes :math:`f_{\mathrm{react}}(c)` of a
reaction mechanism.

In this kind of mechanism, the fluxes are typically defined as a function of the concentration of the components involved in the reaction and the reaction parameters.

In CADET, reaction models are defined for each phase of a unit operation separately, i.e., for the liquid phase and solid phase.
For more information on the interface specification, see :ref:`FFReaction`.

CADET features the following reaction model types:

 - :ref:`mass_action_law_model`
 - :ref:`michaelis_menten_kinetics_model`

It is possible to combine reaction models within the same phase.

Cross phase reaction models
---------------------------

If a reaction is defined between a liquid phase and solid phases, the net fluxes of the reaction are defined
as a function of the concentrations in the respective phases and is called a cross phase reaction.

CADET features the following cross phase reaction types:

 - :ref:`mass_action_law_model_cross_phase`

Further details on the interface specification are provided in :ref:`the corresponding section <FFReaction>`.

Application examples
--------------------

Historically, a chromatography system is modeled as a reaction system without considering any transport phenomenon. We also introduce some reaction-based models that can be solved in CADET:

 - :ref:`thomas_model`
 - :ref:`rate_constant_distribution_theory`

.. _dependence-on-external-function_react:

Dependence on external function
-------------------------------

A reaction model may depend on an external function or profile :math:`T\colon \left[ 0, T_{\mathrm{end}}\right] \times [0, L] \to \mathbb{R}`, where :math:`L` denotes the physical length of the unit operation, or :math:`T\colon \left[0, T_{\mathrm{end}}\right] \to \mathbb{R}` if the unit operation model has no axial length.
By using an external profile, it is possible to account for effects that are not directly modeled in CADET (e.g., temperature).
The dependence of each parameter is modeled by a polynomial of third degree.
For example, the forward rate constant :math:`k_{\mathrm{fwd}}` is really given by

.. math::

    \begin{aligned}
        k_{\mathrm{fwd}}(T) &= k_{\mathrm{fwd},3} T^3 + k_{\mathrm{fwd},2} T^2 + k_{\mathrm{fwd},1} T + k_{\mathrm{fwd},0}.
    \end{aligned}

While :math:`k_{\mathrm{fwd},0}` is set by the original parameter ``XXX_KFWD`` of the file format (``XXX`` being a placeholder for the reaction model), the parameters :math:`k_{\mathrm{fwd},3}`, :math:`k_{\mathrm{fwd},2}`, and :math:`k_{\mathrm{fwd},1}` are given by ``XXX_KFWD_TTT``, ``XXX_KFWD_TT``, and ``XXX_KFWD_T``, respectively.
The identifier of the externally dependent reaction model is constructed from the original identifier by prepending ``EXT_`` (e.g., ``MASS_ACTION_LAW`` is changed into ``EXT_MASS_ACTION_LAW``).
This pattern applies to all parameters and supporting reaction models.
Note that the parameter units have to be adapted to the unit of the external profile by dividing with an appropriate power.

Each parameter of the externally dependent reaction model can depend on a different external source.
The 0-based indices of the external source for each parameter is given in the dataset ``EXTFUN``.
By assigning only one index to ``EXTFUN``, all parameters use the same source.
The ordering of the parameters in ``EXTFUN`` is given by the ordering in the file format specification.

For the configuration of external function dependence and more information on model parameters required to define in CADET file format, see section :ref:`FFReaction`.

.. toctree::
    :hidden:
    :glob:

    *
