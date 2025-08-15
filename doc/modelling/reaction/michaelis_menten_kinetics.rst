.. _michaelis_menten_kinetics_model:

Michaelis Menten kinetics

Implements Michaelis-Menten reaction kinetics of the form

.. math::

    \begin{aligned}
        f_\text{react} = S \nu,
    \end{aligned}

where :math:`S` is the stoichiometric matrix and :math:`\nu` a flux vector with

.. math::

    \begin{aligned}
        \nu_{j} = \prod_{i = 1}^{N_{sub}} \nu_{i,j} = \prod_{i = 1}^{N_{sub}} \frac{ v_{\mathrm{max},j} \, c_i}{K_{\mathrm{M},j} + c_i}.
    \end{aligned}

Here,

- :math:`\nu_{i,j}` is the flux of reaction :math:`j` with respect to substrate :math:`i`,
- :math:`N_{sub,j}` is the number of substrates for reaction :math:`j`,
- and :math:`c_i` is the concentration of each substrate.

Each substrate contributes multiplicatively to the overall rate.
The selection of which components act as substrates is controlled via the negative entries of the stoichiometric matrix.

In addition, CADET supports three types of inhibition reactions.
In this case the flux :math:`\nu_{i,j}` can be modified as one of the following:

Competitive inhibition
^^^^^^^^^^^^^^^^^^^^^^

.. math::

    \begin{aligned}
        \nu_{i,j} =  \frac{v_{\mathrm{max},j} \, c_i}{K_{\mathrm{M},j}\,(1 + \sum_{k \in \mathcal{I}_j} \frac{c_{i,k}}{K_{i,j,k}}) + c_i},
    \end{aligned}

where
 - :math:`c_i` is the substrate component and :math:`c_{i,k}` is one inhibitor acting on substrate :math:`c_i`.
 - :math:`K_{I_{i,j,k}}` is the inhibition constant with respect to component :math:`c_{i,k}` in reaction :math:`j`.
        If :math:`K_{I_{i,j,k}} \leq 0`, component :math:`c_{i,k}` does not act as an inhibitor.
 - :math:`\mathcal{I}_j` is the index set of inhibitors in reaction :math:`j`, i.e all the indices :math:`k` where :math:`K_{I_{i,j,k}} > 0`.

Uncompetitive inhibition
^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

    \begin{aligned}
        \nu_{i,j} = \frac{v_{\mathrm{max},j} \, c_i}{K_{\mathrm{M},j} + c_i \, (1 + \sum_{k \in \mathcal{I}_j} \frac{c_{i,k}}{\tilde{K}_{i,j,k}})}
    \end{aligned}

where
 - :math:`c_i` is the substrate component and :math:`c_{i,k}` is one inhibitor acting on substrate :math:`c_i`.
 - :math:`\tilde{K}_{I_{i,j,k}}` is the inhibition constant with respect to component :math:`c_{i,k}` in reaction :math:`j`.
    If :math:`\tilde{K}_{I_{i,j,k}} \leq 0`, component :math:`c_{i,k}` does not act as an inhibitor.
 - :math:`\mathcal{I}_j` is the index set of inhibitors in reaction :math:`j`, i.e all the indices :math:`k` where :math:`\tilde{K}_{I_{i,j,k}} > 0`.

Non-competitive inhibition
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

    \begin{aligned}
       \nu_{i,j} =  \frac{v_{\mathrm{max},j} \, c_i}{(K_{\mathrm{M},j} + c_i) \,(1 + \sum_{k \in \mathcal{I}_j} \frac{c_{i,k}}{K_{i,j,k}})}
    \end{aligned}

where
 - :math:`c_i` is the substrate component and :math:`c_{i,k}` one the inhibitor acting on substrate :math:`c_i`.
 - :math:`K_{I_{i,j,k}}` is the inhibition constant with respect to component :math:`c_{i,k}` in reaction :math:`j`.
        If :math:`K_{I_{i,j,k}} \leq 0`, component :math:`c_{i,k}` does not act as an inhibitor.
 - :math:`\mathcal{I}_j` is the index set of inhibitors in reaction :math:`j`, i.e all the indices :math:`k` where :math:`K_{I_{i,j,k}} > 0`.
Non-competitive inhibition can be viewed as a combination of competitive and uncompetitive inhibition mechanisms, depending on the values of the corresponding inhibition constants.
If the inhibition constant for uncompetitive inhibition :math:`\tilde{K}_{I_{i,j,k}}` is set to the same positive value as the inhibition constant for the
competitive inhibition :math:`K_{I_{i,j,k}}`, the inhibition type is non-competitive.

For configuration information please refer to :ref:`michaelis_menten_kinetics_config`.


Monod kinetics
^^^^^^^^^^^^^^

The Michaelis-Menten model in CADET can also be used to represent Monod kinetics, as the mathematical formulation is identical.
The Monod equation is commonly used to describe microbial growth processes and corresponds to Michaelis-Menten kinetics with only one substrate.
By choosing appropriate parameters, Monod kinetics can be directly simulated with this model.


Literature
^^^^^^^^^^
- Segel, I. H. (1993). Enzyme kinetics: Behavior and analysis of rapid equilibrium and steady-state enzyme systems. John Wiley & Sons.
