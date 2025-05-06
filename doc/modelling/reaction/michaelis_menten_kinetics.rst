.. _michaelis_menten_kinetics_model_with_inhibition:

Michaelis Menten kinetics with inhibition
-------------------------

Implements liquid phase Michaelis-Menten reaction kinetics of the form

.. math::
    \begin{aligned}
        f_\text{react} = S \mathbf{\nu},
    \end{aligned}

where :math:`S` is the stoichiometric matrix and the fluxes are given by

.. math::

    \begin{aligned}
        \nu_j = \prod_i^{N_{sub}} \frac{\mu_{\mathrm{max},j} \, c_i}{k_{\mathrm{MM},j} + c_i},
    \end{aligned}

where

- :math:`j` is the reaction index,
- :math:`N_{sub}` is the number of substrates in the reaction,
- :math:`c_i` is the substrate component,
- :math:`\mu_{\mathrm{max},j}`, is the limiting rate approached by the system at saturation,
- :math:`k_{\mathrm{MM},j}` is the Michaelis constant, which is defined as the concentration of substrate at which the reaction rate is half of :math:`\mu_{\mathrm{max},j}`.


In addition, a substrate :math:`i` might be inhibited by other components.
There are three types of inhibition implemented.

Competitive inhibition
^^^^^^^^^^^^^^^^^^^^^^^

.. math::

    \begin{aligned}
        \nu_{i,j} =  \frac{\mu_{\mathrm{max},j} \cdot c_i}{k_{\mathrm{MM},j}(1+\sum_k^{N_{inh}} \frac{c_{i,k}}{K_{i,j,k}}) + c_i},
    \end{aligned}

where
 - :math:`c_i` is the substrate component and :math:`c_{i,k}` is the inhibitor  acting on substrate :math:`c_i`.
 - :math:`N_{inh}` is the number of inhibitors in the reaction,
 - :math:`K_{i,j,k}` is the inhibition constant with respect to component :math:`c_{i,k}` in reaction :math:`j`.
        If :math:`K_{i,j,k} \leq 0`, component :math:`c_{i,k}` does not act as an inhibitor.


Uncompetitive inhibition
^^^^^^^^^^^^^^^^^^^^^^^

.. math::

    \begin{aligned}
        \nu_{i,j} = \frac{\mu_{\mathrm{max},j} \, c_i}{k_{\mathrm{MM},j} + c_i \, (1 + \sum_k^{N_{inh}} \frac{c_{i,k}}{\tilde{K}_{i,j,k}})}
    \end{aligned}

where
 - :math:`c_i` is the substrate component and :math:`c_{i,k}` is the inhibitor  acting on substrate :math:`c_i`.
 - :math:`N_{inh}` is the number of inhibitors in the reaction,
 - :math:`\tilde{K}_{i,j,k}` is the inhibition constant with respect to component :math:`c_{i,k}` in reaction :math:`j`.
    If :math:`\tilde{K}_{i,j,k} \leq 0`, component :math:`c_{i,k}` does not act as an inhibitor.

Non-competitive inhibition
^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

    \begin{aligned}
       \nu_{i,j} =  \frac{\mu_{\mathrm{max},j} \, c_i}{(k_{\mathrm{MM},j} + c_i) \,(1 + \sum_k^{N_{inh}} \frac{c_{i,k}}{K_{i,j,k}})}
    \end{aligned}

where
 - :math:`c_i` is the substrate component and :math:`c_{i,k}` is the inhibitor acting on substrate :math:`c_i`.
 - :math:`N_{inh}` is the number of inhibitors in the reaction,
 - :math:`K_{i,j,k}` is the inhibition constant with respect to component :math:`c_{i,k}` in reaction :math:`j`.
        If :math:`K_{i,j,k} \leq 0`, component :math:`c_{i,k}` does not act as an inhibitor.
Non-competitive inhibition can be viewed as a combination of competitive and uncompetitive inhibition mechanisms, depending on the values of the corresponding inhibition constants.
If the inhibition constant for uncompetitive inhibition :math:`\tilde{K}_{i,j,k}` is set to the same positive value as the inhibition constant for the 
competitive inhibition :math:`K_{i,j,k}`, the inhibition type is non-competitive.

The rate for reaction :math:`j` is given by the product of the individual reaction rates of each substrate
    
.. math::

    \begin{aligned}
        \nu_j = \prod_{i=1}^{N_{sub}} \nu_{i,j}.
    \end{aligned}
