.. _michaelis_menten_kinetics_model:

Michaelis Menten kinetics
-------------------------

Implements liquid phase Michaelis-Menten reaction kinetics of the form

.. math::

    \begin{aligned}
        f_\text{react} = S \mathbf{\nu},
    \end{aligned}

where :math:`S` is the stoichiometric matrix and the fluxes are given by

.. math::

    \begin{aligned}
        \nu_j = \mu_{\mathrm{max},j} \prod_{i \in \mathcal{S}_j} \frac{{c_i}}{k_{\mathrm{MM},j} + c_i},
    \end{aligned}

where

- :math:`j` is the reaction index,
- :math:`\mathcal{S}_j` is the set of substrate indices for reaction :math:`j`,
- :math:`\mu_{\mathrm{max},j}`, is the limiting rate approached by the system at saturation,
- :math:`k_{\mathrm{MM},j}` is the Michaelis constant, which is defined as the concentration of substrate at which the reaction rate is half of :math:`\mu_{\mathrm{max},j}`.


In addition, the reaction might be inhibited by other components due to noncompetitive inhibition.
In this case, the flux has the form

.. math::

    \begin{aligned}
        \nu_j = \mu_{\mathrm{max},j} \prod_{i \in \mathcal{S}_j}\frac{c_i}{k_{\mathrm{MM},j} + c_i} \cdot \frac{1}{1 + \sum_{k = 1}^{n} \frac{c_{k}}{/k_{I,j,k}}}.
    \end{aligned}

Here, :math: `k_{\mathrm{I},j,k}` is the inhibition constant with respect to component :math:`k` in reaction :math:`j`.
If :math: `k_{\mathrm{I},j,k} \leq 0`, component :math:`k` does not act as an inhibitor.

Note: Currently, the model does not allow substrates to function as inhibitors.

References
^^^^^^^^^^
.. [1] Irwin H. Segel (1975). Enzyme Kinetics: Behavior and Analysis of Rapid Equilibrium and Steady State Enzyme Systems.
