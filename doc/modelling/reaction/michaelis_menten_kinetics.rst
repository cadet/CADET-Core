.. _michaelis_menten_kinetics_model:

Michaelis Menten kinetics
=========================

Implements Michaelis-Menten reaction kinetics of the form

.. math::

    \begin{aligned}
        f_\text{react} = S \nu,
    \end{aligned}

where :math:`S` is the stoichiometric matrix and :math:`\nu` a flux vector with

.. math::

    \begin{aligned}
        \nu_{j} = v_{\mathrm{max},j} \prod_{i = 1}^{N_{sub,j}} \nu_{i,j} = v_{\mathrm{max},j} \prod_{i = 1}^{N_{sub,j}} \frac{ c_{i,j}}{K_{\mathrm{M},i,j} + c_{i,j}}
    \end{aligned}

where

- :math:`\nu_{j}` is the flux of reaction :math:`j` and :math:`\nu_{i,j}` is the flux of reaction :math:`j` with respect to substrate :math:`i`,
- :math:`N_{sub,j}` is the number of substrates for reaction :math:`j`,
- :math:`v_{\mathrm{max},j}` is the maximum reaction rate in reaction :math:`j`,
- and :math:`c_{i,j}` is the concentration of substrate :math:`i` in reaction :math:`j`.
- :math:`K_{\mathrm{M}_{i,j}}` is the Michaelis constant for substrate :math:`i` in reaction :math:`j`.

The selection of which components act as substrates is controlled via the sign of the entries of the stoichiometric matrix.
For more information about the configuration can be found in :ref:`michaelis_menten_kinetics_config`.

In addition, CADET supports three types of inhibition reactions.
In this case the flux :math:`\nu_{i,j}` can be modified as one of the following:

Competitive Inhibition
^^^^^^^^^^^^^^^^^^^^^^^
In competitive inhibition, the inhibitor binds at the enzyme's active site. The modified flux expression is:

.. math::

    \begin{aligned}
        \nu_{i,j} =  \frac{ c_{i,j}}{K_{\mathrm{M},i,j}\,(1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{K^{c}_{I_{k}}}) + c_{i,j}},
    \end{aligned}

where
 - :math:`c_{i,j}` is the substrate component and :math:`c_{k}` is one inhibitor acting on substrate :math:`c_{i,j}`,
 - :math:`K^{c}_{I_{k}}` is the inhibition constant with respect to inhibitor :math:`c_{k}` i.e if :math:`K^{c}_{I_{k}} > 0`, component :math:`c_{k}` acts as an inhibitor to substrate :math:`c_{i,j}`,
 - :math:`\mathcal{I}_{i,j}` is the index set of inhibitors for substrate :math:`c_{i,j}`, i.e the indices :math:`k` where :math:`K^{c}_{I_{k}} > 0`.

Uncompetitive Inhibition
^^^^^^^^^^^^^^^^^^^^^^^^

In an uncompetitive inhibition, the inhibitor binds to the enzyme-substrate complex, preventing the reaction from proceeding. The modified flux expression is:

.. math::

    \begin{aligned}
        \nu_{i,j} = \frac{c_{i,j}}{K_{\mathrm{M},i,j} + c_{i,j} \, (1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{\tilde{K}_{k}})},
    \end{aligned}

where
 - :math:`c_{i,j}` is the substrate component and :math:`c_{k}` is one inhibitor acting on substrate :math:`c_{i,j}`.
 - :math:`K^{u}_{I_{k}}` is the inhibition constant with respect to component :math:`c_{k}` in reaction :math:`j` i.e if :math:`K^{u}_{I_{k}} > 0`, component :math:`c_{k}` acts as an inhibitor to substrate :math:`c_{i,j}`.
 - :math:`\mathcal{I}_{i,j}` is the index set of inhibitors in reaction :math:`j`, i.e the indices :math:`k` where :math:`K^{u}_{I_{k}} > 0`.

Mixed Inhibition
^^^^^^^^^^^^^^^^

In non-competitive inhibition, the inhibitor can bind to both the enzyme and the enzyme-substrate complex, preventing the reaction from proceeding. The modified flux expression is:

.. math::

    \begin{aligned}
       \nu_{i,j} =  \frac{c_{i,j}}{ K_{\mathrm{M},i,j} \,(1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{K^{c}_{I_{k}}}) + c_{i,j} \,(1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{K^{u}_{I_{k}}})},
    \end{aligned}

where
 - :math:`c_{i,j}` is the substrate component and :math:`c_{k}` is one the inhibitor acting on substrate :math:`c_{i,j}`.
 - :math:`K^{c}_{I_{k}}`and :math:`K^{u}_{I_{k}}` are the inhibition constants with respect to component :math:`c_{k}` in reaction :math:`j` i.e if :math:`K^{c}_{I_{k}} > 0`, component :math:`c_{k}` acts as an inhibitor to substrate :math:`c_{i,j}`.
 - :math:`\mathcal{I}_{i,j}` is the index set of inhibitors in reaction :math:`j`, i.e the indices :math:`k` where :math:`K^{c}_{I_{k}} > 0`.


Non-Competitive Inhibition
^^^^^^^^^^^^^^^^^^^^^^^^^^

Non-competitive inhibition is a form of mixed inhibition where the inhibitor binds to both the enzyme and the enzyme-substrate complex with the **same affinity**.

.. math::

    \begin{aligned}
       \nu_{i,j} =  \frac{c_{i,j}}{(K_{\mathrm{M},i,j} + c_{i,j}) \,(1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{K^{n}_{I_{k}}})}
    \end{aligned}

where
 - :math:`c_{i,j}` is the substrate component and :math:`c_{k}` is one the inhibitor acting on substrate :math:`c_{i,j}`.
 - :math:`K^{n}_{I_{k}}` is the inhibition constant with respect to component :math:`c_{k}` in reaction :math:`j` i.e if :math:`K^{n}_{I_{k}} > 0`, component :math:`c_{k}` acts as an inhibitor to substrate :math:`c_{i,j}`.
 - :math:`\mathcal{I}_{i,j}` is the index set of inhibitors in reaction :math:`j`, i.e the indices :math:`k` where :math:`K^{n}_{I_{k}} > 0`.
Note that the inhibition constant for the non-competitive inhibition is indirectly given if :math:` K^{c}_{I_{k}} = K^{u}_{I_{k}} = K^{n}_{I_{k}}`

For configuration information please refer to :ref:`michaelis_menten_kinetics_config`.


Monod Kinetics
^^^^^^^^^^^^^^

The Michaelis-Menten model in CADET can also be used to represent Monod kinetics, as the mathematical formulation is identical.
The Monod equation is commonly used to describe microbial growth processes and corresponds to Michaelis-Menten kinetics with only one substrate.

The Monod equation for microbial growth has the form:

.. math::

    \begin{aligned}
        \mu = \frac{\mu_{\mathrm{max}} \, c_S}{K_S + c_S}
    \end{aligned}

where:

- :math:`\mu` is the specific growth rate
- :math:`\mu_{\mathrm{max}}` is the maximum specific growth rate
- :math:`c_S` is the substrate concentration
- :math:`K_S` is the saturation constant (half-saturation constant)

By choosing a Michaelis-Menten kinetics configuration with one substrate and setting the the parameter accordingly,
i.e :math:`\mu_{\mathrm{max}} = v_{\mathrm{max}}`, :math:`K_S = K_{\mathrm{M},0,0}` and :math:`c_S = c_{0,0}`,
the Monod equation can be expressed in the same form as the Michaelis-Menten kinetics.

Literature
^^^^^^^^^^

- Segel, I. H. (1993). Enzyme kinetics: Behavior and analysis of rapid equilibrium and steady-state enzyme systems. John Wiley & Sons.
- Monod, Jacques. 1949. “The Growth of Bacterial Cultures.” *Annual Review of Microbiology* 3 (1): 371–394. https://doi.org/10.1146/annurev.mi.03.100149.002103
