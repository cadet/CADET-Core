.. _mass_action_law_model:

Mass action law
---------------

The mass action law reaction model is suitable for most reactions.
Note that the concentrations are directly used for calculating the fluxes.
Hence, the model only holds for dilute solutions under the assumption of a well-stirred reaction vessel.
These assumptions can be weakened by passing to the generalized mass action law, which uses chemical activities instead of concentrations.

The mass action law states that the speed of a reaction is proportional to the product of the concentrations of their reactants.
The net flux for component :math:`i` is given by

.. math::

    \begin{aligned}
        f_{\mathrm{react},i}^l\left(c^l\right) &= \sum_{j=0}^{N_{\mathrm{react}}-1} s_{i,j}^l \varphi^l_j\left(c^l\right), \\
        \varphi^l_j(c^l) &= k^l_{\mathrm{fwd},j} \prod_{\ell=0}^{N_{\mathrm{comp}}-1} \left(c^l_{\ell}\right)^{e^l_{\mathrm{fwd},\ell,j}} - k^l_{\mathrm{bwd},j} \prod_{\ell=0}^{N_{\mathrm{comp}}-1} \left(c^l_{\ell}\right)^{e^l_{\mathrm{bwd},\ell,j}},
    \end{aligned}

where :math:`S^l = (s^l_{i,j}) \in \mathbb{R}^{N_{\mathrm{comp}} \times N_{\mathrm{react}}}` is the stoichiometric matrix, :math:`\varphi^l_j(c)` is the net flux of reaction :math:`j`, and :math:`k^l_{\mathrm{fwd},j}` and :math:`k^l_{\mathrm{bwd},j}` are the rate constants.
The matrices :math:`E^l_{\mathrm{fwd}} = (e^l_{\mathrm{fwd},\ell,j}) \in \mathbb{R}^{N_{\mathrm{comp}} \times N_{\mathrm{react}}}` and :math:`E^l_{\mathrm{bwd}} = (e^l_{\mathrm{bwd},\ell,j}) \in \mathbb{R}^{N_{\mathrm{comp}} \times N_{\mathrm{react}}}` are usually derived by the order of the reaction, that is,

.. math::
    :label: MRMassActionLawExpMatDefault

    \begin{aligned}
        e^l_{\mathrm{fwd},\ell,j} &= \max(0, -s^l_{\ell,j}), \\
        e^l_{\mathrm{bwd},\ell,j} &= \max(0, s^l_{\ell,j}).
    \end{aligned}

However, these defaults can be changed by providing those matrices.

In situations where both liquid and solid phase are present (e.g., in a bead), the respective other phase may act as a modifier in the net flux equation.
For example, consider reactions in the liquid phase of a particle given by

.. math::

    \begin{aligned}
        f_{\mathrm{react},i}^p\left(c^p, c^s\right) &= \sum_{j=0}^{N_{\mathrm{react}}-1} s_{i,j}^p \varphi^p_j\left(c^p, c^s\right),\end{aligned}

where

.. math::

    \begin{split}
        \varphi^p_j(c^p, c^s) = k^p_{\mathrm{fwd},j} &\left[\prod_{\ell=0}^{N_{\mathrm{comp}}-1} \left(c^p_{\ell}\right)^{e^p_{\mathrm{fwd},\ell,j}}\right] \left[\prod_{m=0}^{\sum_{i=0}^{N_{\mathrm{comp}}-1} N_{\mathrm{bnd},i}-1} \left(c^s_{m}\right)^{e^{ps}_{\mathrm{fwd},m,j}}\right] \\
         - k^p_{\mathrm{bwd},j} &\left[\prod_{\ell=0}^{N_{\mathrm{comp}}-1} \left(c^p_{\ell}\right)^{e^p_{\mathrm{bwd},\ell,j}}\right] \left[\prod_{m=0}^{\sum_{i=0}^{N_{\mathrm{comp}}-1} N_{\mathrm{bnd},i}-1} \left(c^s_{m}\right)^{e^{ps}_{\mathrm{bwd},m,j}}\right].
    \end{split}

The forward and backward rates of the liquid phase particle reactions can be modified by a power of every bound state in the solid phase of the particle.
The exponents of these powers are given by the matrices :math:`E^{ps}_{\mathrm{fwd}} = (e^{ps}_{\mathrm{fwd},m,j})` and :math:`E^{ps}_{\mathrm{bwd}} = (e^{ps}_{\mathrm{bwd},m,j})`, which are both of size :math:`(\sum_i N_{\mathrm{bnd},i}) \times N_{\mathrm{react}}`.
Whereas the exponent matrices :math:`E^{p}_{\mathrm{fwd}}, E^{p}_{\mathrm{bwd}} \in \mathbb{R}^{N_{\mathrm{comp}} \times N_{\mathrm{react}}}` are initialized based on the stoichiometric matrix :math:`S^{p} \in \mathbb{R}^{N_{\mathrm{comp}} \times N_{\mathrm{react}}}`, see Eq. :eq:`MRMassActionLawExpMatDefault`, the exponent matrices :math:`E^{ps}_{\mathrm{fwd}}, E^{ps}_{\mathrm{bwd}}` of the modifier terms default to :math:`0`.

Vice versa, the rates of solid phase reactions can be modified by liquid phase concentrations.
The corresponding exponent matrices :math:`E^{sp}_{\mathrm{fwd}} = (e^{sp}_{\mathrm{fwd},\ell,j})` and :math:`E^{sp}_{\mathrm{bwd}} = (e^{sp}_{\mathrm{bwd},\ell,j})` are both of size :math:`N_{\mathrm{comp}} \times N_{\mathrm{react}}`.

.. math::

    \begin{aligned}
        f_{\mathrm{react},i}^s\left(c^s, c^p\right) &= \sum_{j=0}^{N_{\mathrm{react}}-1} s_{i,j}^s \varphi^s_j\left(c^s, c^p\right),
    \end{aligned}

where

.. math::

    \begin{split}
        \varphi^s_j(c^s, c^p) = k^s_{\mathrm{fwd},j} &\left[\prod_{m=0}^{\sum_{i=0}^{N_{\mathrm{comp}}-1} N_{\mathrm{bnd},i}-1} \left(c^s_{m}\right)^{e^{s}_{\mathrm{fwd},m,j}}\right] \left[\prod_{\ell=0}^{N_{\mathrm{comp}}-1} \left(c^p_{\ell}\right)^{e^{sp}_{\mathrm{fwd},\ell,j}}\right] \\
        - k^p_{\mathrm{bwd},j} &\left[\prod_{m=0}^{\sum_{i=0}^{N_{\mathrm{comp}}-1} N_{\mathrm{bnd},i}-1} \left(c^s_{m}\right)^{e^{s}_{\mathrm{bwd},m,j}}\right] \left[\prod_{\ell=0}^{N_{\mathrm{comp}}-1} \left(c^p_{\ell}\right)^{e^{sp}_{\mathrm{bwd},\ell,j}}\right].
    \end{split}

Whereas the exponent matrices :math:`E^{s}_{\mathrm{fwd}}, E^{s}_{\mathrm{bwd}} \in \mathbb{R}^{(\sum_i N_{\mathrm{bnd},i}) \times N_{\mathrm{react}}}` are initialized based on the stoichiometric matrix :math:`S^{s} \in \mathbb{R}^{(\sum_i N_{\mathrm{bnd},i}) \times N_{\mathrm{react}}}`, see Eq. :eq:`MRMassActionLawExpMatDefault`, the exponent matrices :math:`E^{sp}_{\mathrm{fwd}}, E^{sp}_{\mathrm{bwd}}` of the modifier terms default to :math:`0`.


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

For more information on model parameters required to define in CADET file format, see :ref:`mass_action_law_config`.
