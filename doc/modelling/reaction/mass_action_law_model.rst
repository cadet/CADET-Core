.. _mass_action_law_model:

Mass Action Law
---------------

The mass action law reaction model is suitable for most reactions.
Note that the concentrations are directly used for calculating the fluxes.
Hence, the model only holds for dilute solutions under the assumption of a well-stirred reaction vessel.
These assumptions can be weakened by passing to the generalized mass action law, which uses chemical activities instead of concentrations.

The mass action law states that the speed of a reaction is proportional to the product of the concentrations of their reactants.

The indices used in the equations have the following meaning:

- :math:`i` is the component index, ranging from :math:`0` to :math:`N_{\mathrm{comp}}-1`
- :math:`j` is the reaction index, ranging from :math:`0` to :math:`N_{\mathrm{react}}-1`
- :math:`\ell` is the component index used in products and sums, ranging from :math:`0` to :math:`N_{\mathrm{comp}}-1`

The net flux for component :math:`i` is given by

.. math::

    \begin{aligned}
        f_{\mathrm{react},i}\left(c\right) &= \sum_{j=0}^{N_{\mathrm{react}}-1} s_{i,j} \varphi_j\left(c\right), \\
        \varphi_j(c) &= k_{\mathrm{fwd},j} \prod_{\ell=0}^{N_{\mathrm{comp}}-1} \left(c_{\ell}\right)^{e_{\mathrm{fwd},\ell,j}} - k_{\mathrm{bwd},j} \prod_{\ell=0}^{N_{\mathrm{comp}}-1} \left(c_{\ell}\right)^{e_{\mathrm{bwd},\ell,j}},
    \end{aligned}

where :math:`S = (s_{i,j}) \in \mathbb{R}^{N_{\mathrm{comp}} \times N_{\mathrm{react}}}` is the stoichiometric matrix, :math:`\varphi_j(c)` is the net flux of reaction :math:`j`, and :math:`k_{\mathrm{fwd},j}` and :math:`k_{\mathrm{bwd},j}` are the rate constants.
The matrices :math:`E_{\mathrm{fwd}} = (e_{\mathrm{fwd},\ell,j}) \in \mathbb{R}^{N_{\mathrm{comp}} \times N_{\mathrm{react}}}` and :math:`E_{\mathrm{bwd}} = (e_{\mathrm{bwd},\ell,j}) \in \mathbb{R}^{N_{\mathrm{comp}} \times N_{\mathrm{react}}}` are usually derived by the order of the reaction, that is,

.. math::
    :label: MRMassActionLawExpMatDefault

    \begin{aligned}
        e_{\mathrm{fwd},\ell,j} &= \max(0, -s_{\ell,j}), \\
        e_{\mathrm{bwd},\ell,j} &= \max(0, s_{\ell,j}).
    \end{aligned}

However, these defaults can be changed by providing those matrices.

This model applies to reactions occurring within a single phase (liquid, particle liquid, or solid phase) only.

For reactions involving multiple phases, see :ref:`mass_action_law_model_cross_phase`.


Correlation of forward- and backward rate constants
---------------------------------------------------

Note that forward rate constant :math:`k_{\mathrm{fwd},i}` and backward
rate constant :math:`k_{\mathrm{bwd},i}` of reaction :math:`i` are
linearly correlated due to the form of the equilibrium constant
:math:`k_{\mathrm{eq},i}`:

.. math::

    \begin{aligned}
        k_{\mathrm{fwd},i} = k_{\mathrm{eq},i} k_{\mathrm{bwd},i}.
    \end{aligned}

This correlation can potentially degrade performance of some optimization algorithms.
The parameters can be decoupled by reparameterization:

.. math::

    \begin{aligned}
        r_{\mathrm{net},i} &= k_{\mathrm{fwd},i} f_{\mathrm{fwd},i} - k_{\mathrm{bwd},i} f_{\mathrm{bwd},i}\\
        &= k_{\mathrm{bwd},i} \left[ k_{\mathrm{eq},i} f_{\mathrm{fwd},i} - f_{\mathrm{bwd},i} \right] \\
        &= k_{\mathrm{fwd},i} \left[ f_{\mathrm{fwd},i} - \frac{1}{k_{\mathrm{eq},i}} f_{\mathrm{bwd},i} \right].
    \end{aligned}

This can be achieved by a (nonlinear) parameter transform

.. math::

    \begin{aligned}
        F\left( k_{\mathrm{eq},i}, k_{\mathrm{bwd},i} \right) &= \begin{pmatrix} k_{\mathrm{eq},i} k_{\mathrm{bwd},i} \\ k_{\mathrm{bwd},i} \end{pmatrix} \\
        \text{ with Jacobian } J_F\left( k_{\mathrm{eq},i}, k_{\mathrm{bwd},i} \right) &= \begin{pmatrix} k_{\mathrm{bwd},i} & k_{\mathrm{eq},i} \\ 0 & 1 \end{pmatrix}.
    \end{aligned}

For more information on model parameters required to define in CADET file format, see :ref:`_mass_action_law_config`.
