.. _mass_action_law_model_cross_phase:

Mass Action Law Cross Phase
---------------------------

The mass action law cross phase reaction model describes reactions that occur between the liquid and solid phases.

This model extends the basic mass action law to handle interactions between phases, such as reactions between liquid phase components and bound states in the solid phase.

As for the singe phase Mass Action Law note, that the concentrations are directly used for calculating the fluxes.
Hence, the model only holds for dilute solutions under the assumption of a well-stirred reaction vessel.
These assumptions can be weakened by passing to the generalized mass action law, which uses chemical activities instead of concentrations.

The mass action law states that the speed of a reaction is proportional to the product of the concentrations of their reactants.

For single-phase reactions within a single phase only, see :ref:`mass_action_law_model`.

Cross-Phase Reactions
~~~~~~~~~~~~~~~~~~~~~

In situations where both liquid and solid phase are present (e.g., in a bead), reactions can occur between phases, where one phase may act as a modifier in the net flux equation of the other phase.

In the following :math:`l` denotes the liquid phase and :math:`s` the solid phase.
Depending on the unit operation, the liquid phase can be considered as the bulk or particle phase.

The indices used in the equations have the following meaning:

- :math:`i` is the component index,
- :math:`j` is the reaction index, ranging from :math:`0` to :math:`N_{\mathrm{react}}-1`
- :math:`\ell` is the component index for liquid phase components, ranging from :math:`0` to :math:`N_{\mathrm{comp}}-1`
- :math:`m` is the bound state index for solid phase components, ranging from :math:`0` to :math:`\sum_{i=0}^{N_{\mathrm{comp}}-1} N_{\mathrm{bnd},i}-1`

Liquid Phase Reactions Modified by Solid Phase
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The exponent matrices are usually derived by the order of the reaction, that is,

.. math::
    :label: MRMassActionLawExpMatDefaultCrossPhase

    \begin{aligned}
        e^l_{\mathrm{fwd},\ell,j} &= \max(0, -s^l_{\ell,j}), \\
        e^l_{\mathrm{bwd},\ell,j} &= \max(0, s^l_{\ell,j}).
    \end{aligned}

However, these defaults can be changed by providing those matrices.

For example, consider reactions in the liquid phase of a particle given by

.. math::

    \begin{aligned}
        f_{\mathrm{react},i}^l\left(c^l, c^s\right) &= \sum_{j=0}^{N_{\mathrm{react}}-1} s_{i,j}^l \varphi^l_j\left(c^l, c^s\right),\end{aligned}

where

.. math::

    \begin{split}
        \varphi^l_j(c^l, c^s) = k^l_{\mathrm{fwd},j} &\left[\prod_{\ell=0}^{N_{\mathrm{comp}}-1} \left(c^l_{\ell}\right)^{e^l_{\mathrm{fwd},\ell,j}}\right] \left[\prod_{m=0}^{\sum_{i=0}^{N_{\mathrm{comp}}-1} N_{\mathrm{bnd},i}-1} \left(c^s_{m}\right)^{e^{ls}_{\mathrm{fwd},m,j}}\right] \\
         - k^l_{\mathrm{bwd},j} &\left[\prod_{\ell=0}^{N_{\mathrm{comp}}-1} \left(c^l_{\ell}\right)^{e^l_{\mathrm{bwd},\ell,j}}\right] \left[\prod_{m=0}^{\sum_{i=0}^{N_{\mathrm{comp}}-1} N_{\mathrm{bnd},i}-1} \left(c^s_{m}\right)^{e^{ls}_{\mathrm{bwd},m,j}}\right].
    \end{split}

The forward and backward rates of the liquid phase particle reactions can be modified by a power of every bound state in the solid phase of the particle.
The exponents of these powers are given by the matrices :math:`E^{ls}_{\mathrm{fwd}} = (e^{ls}_{\mathrm{fwd},m,j})` and :math:`E^{ls}_{\mathrm{bwd}} = (e^{ls}_{\mathrm{bwd},m,j})`, which are both of size :math:`(\sum_i N_{\mathrm{bnd},i}) \times N_{\mathrm{react}}`.
Whereas the exponent matrices :math:`E^{p}_{\mathrm{fwd}}, E^{p}_{\mathrm{bwd}} \in \mathbb{R}^{N_{\mathrm{comp}} \times N_{\mathrm{react}}}` are initialized based on the stoichiometric matrix :math:`S^{p} \in \mathbb{R}^{N_{\mathrm{comp}} \times N_{\mathrm{react}}}`, see Eq. :eq:`MRMassActionLawExpMatDefaultCrossPhase`, the exponent matrices :math:`E^{ls}_{\mathrm{fwd}}, E^{ls}_{\mathrm{bwd}}` of the modifier terms default to :math:`0`.

Solid Phase Reactions Modified by Liquid Phase
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Vice versa, the rates of solid phase reactions can be modified by liquid phase concentrations.
The corresponding exponent matrices :math:`E^{sp}_{\mathrm{fwd}} = (e^{sp}_{\mathrm{fwd},\ell,j})` and :math:`E^{sp}_{\mathrm{bwd}} = (e^{sp}_{\mathrm{bwd},\ell,j})` are both of size :math:`N_{\mathrm{comp}} \times N_{\mathrm{react}}`.

.. math::

    \begin{aligned}
        f_{\mathrm{react},i}^s\left(c^s, c^l\right) &= \sum_{j=0}^{N_{\mathrm{react}}-1} s_{i,j}^s \varphi^s_j\left(c^s, c^l\right),
    \end{aligned}

where

.. math::

    \begin{split}
        \varphi^s_j(c^s, c^l) = k^s_{\mathrm{fwd},j} &\left[\prod_{m=0}^{\sum_{i=0}^{N_{\mathrm{comp}}-1} N_{\mathrm{bnd},i}-1} \left(c^s_{m}\right)^{e^{s}_{\mathrm{fwd},m,j}}\right] \left[\prod_{\ell=0}^{N_{\mathrm{comp}}-1} \left(c^l_{\ell}\right)^{e^{sp}_{\mathrm{fwd},\ell,j}}\right] \\
        - k^l_{\mathrm{bwd},j} &\left[\prod_{m=0}^{\sum_{i=0}^{N_{\mathrm{comp}}-1} N_{\mathrm{bnd},i}-1} \left(c^s_{m}\right)^{e^{s}_{\mathrm{bwd},m,j}}\right] \left[\prod_{\ell=0}^{N_{\mathrm{comp}}-1} \left(c^l_{\ell}\right)^{e^{sp}_{\mathrm{bwd},\ell,j}}\right].
    \end{split}

Whereas the exponent matrices :math:`E^{s}_{\mathrm{fwd}}, E^{s}_{\mathrm{bwd}} \in \mathbb{R}^{(\sum_i N_{\mathrm{bnd},i}) \times N_{\mathrm{react}}}` are initialized based on the stoichiometric matrix :math:`S^{s} \in \mathbb{R}^{(\sum_i N_{\mathrm{bnd},i}) \times N_{\mathrm{react}}}`, see Eq. :eq:`MRMassActionLawExpMatDefaultCrossPhase`, the exponent matrices :math:`E^{sp}_{\mathrm{fwd}}, E^{sp}_{\mathrm{bwd}}` of the modifier terms default to :math:`0`.

For more information on model parameters required to define in CADET file format, see :ref:`mass_action_law_cross_phase_config`.
