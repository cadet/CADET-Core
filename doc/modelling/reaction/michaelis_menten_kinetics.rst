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
        \nu_{j} = v_{\mathrm{max},j} \prod_{i = 1}^{N_{sub,j}} \nu_{i,j} = v_{\mathrm{max},j} \prod_{i = 1}^{N_{sub,j}} \frac{ c_{i,j}}{K_{\mathrm{M}_{i,j}} + c_{i,j}}
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
        \nu_{i,j} =  \frac{c_{i,j}}{K_{\mathrm{M}_{i,j}}\,(1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{K^{c}_{I_{k}}}) + c_{i,j}},
    \end{aligned}

where
 - :math:`c_{i,j}` is the substrate component and :math:`c_{k}` is one inhibitor acting on substrate :math:`c_{i,j}`,
 - :math:`K^{c}_{I_{k}}` is the inhibition constant with respect to inhibitor :math:`c_{k}` i.e if :math:`K^{c}_{I_{k}} > 0`, component :math:`c_{k}` acts as an inhibitor to substrate :math:`c_{i,j}`,
 - :math:`\mathcal{I}^{c}_{i,j}` is the index set of inhibitors for substrate :math:`c_{i,j}`, i.e the indices :math:`k` where :math:`K^{c}_{I_{k}} > 0`.

Uncompetitive Inhibition
^^^^^^^^^^^^^^^^^^^^^^^^

In an uncompetitive inhibition, the inhibitor binds to the enzyme-substrate complex, preventing the reaction from proceeding. The modified flux expression is:

.. math::

    \begin{aligned}
        \nu_{i,j} = \frac{c_{i,j}}{K_{\mathrm{M}_{i,j}} + c_{i,j} \, (1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{K^{uc}_{I_{k}}})},
    \end{aligned}

where
 - :math:`c_{i,j}` is the substrate component and :math:`c_{k}` is one inhibitor acting on substrate :math:`c_{i,j}`.
 - :math:`K^{uc}_{I_{k}}` is the inhibition constant with respect to component :math:`c_{k}` in reaction :math:`j` i.e if :math:`K^{uc}_{I_{k}} > 0`, component :math:`c_{k}` acts as an inhibitor to substrate :math:`c_{i,j}`.
 - :math:`\mathcal{I}^{uc}_{i,j}` is the index set of inhibitors in reaction :math:`j`, i.e the indices :math:`k` where :math:`K^{uc}_{I_{k}} > 0`.

Mixed Inhibition
^^^^^^^^^^^^^^^^

In mixed inhibition, the inhibitor can bind to both the enzyme and the enzyme-substrate complex, preventing the reaction from proceeding. The modified flux expression is:

.. math::

    \begin{aligned}
       \nu_{i,j} =  \frac{c_{i,j}}{ K_{\mathrm{M}_{i,j}} \,(1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{K^{c}_{I_{k}}}) + c_{i,j} \,(1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{K^{uc}_{I_{k}}})},
    \end{aligned}

where
 - :math:`c_{i,j}` is the substrate component and :math:`c_{k}` is one the inhibitor acting on substrate :math:`c_{i,j}`.
 - :math:`K^{c}_{I_{k}}` and :math:`K^{uc}_{I_{k}}` are the inhibition constants with respect to component :math:`c_{k}` in reaction :math:`j` i.e if :math:`K^{c}_{I_{k}} > 0`, component :math:`c_{k}` acts as an inhibitor to substrate :math:`c_{i,j}`.
 - :math:`\mathcal{I}^{c}_{i,j}` is the index set of inhibitors in reaction :math:`j`, i.e the indices :math:`k` where :math:`K^{c}_{I_{k}} > 0`,
 - :math:`\mathcal{I}^{uc}_{i,j}` is the index set of inhibitors in reaction :math:`j`, i.e the indices :math:`k` where :math:`K^{uc}_{I_{k}} > 0`.


Non-Competitive Inhibition
^^^^^^^^^^^^^^^^^^^^^^^^^^

Non-competitive inhibition is a form of mixed inhibition where the inhibitor binds to both the enzyme and the enzyme-substrate complex with the **same affinity**.

.. math::

    \begin{aligned}
       \nu_{i,j} =  \frac{c_{i,j}}{(K_{\mathrm{M}_{i,j}} + c_{i,j}) \,(1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{K^{n}_{I_{k}}})}
    \end{aligned}

where
 - :math:`c_{i,j}` is the substrate component and :math:`c_{k}` is one the inhibitor acting on substrate :math:`c_{i,j}`.
 - :math:`K^{n}_{I_{k}}` is the inhibition constant with respect to component :math:`c_{k}` in reaction :math:`j` i.e if :math:`K^{n}_{I_{k}} > 0`, component :math:`c_{k}` acts as an inhibitor to substrate :math:`c_{i,j}`.
 - :math:`\mathcal{I}_{i,j}` is the index set of inhibitors in reaction :math:`j`, i.e the indices :math:`k` where :math:`K^{n}_{I_{k}} > 0`.
Note that the inhibition constant for the non-competitive inhibition is indirectly given if :math:`K^{c}_{I_{k}} = K^{uc}_{I_{k}} = K^{n}_{I_{k}}`

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
i.e :math:`\mu_{\mathrm{max}} = v_{\mathrm{max}}`, :math:`K_S = K_{\mathrm{M}_{0,0}}` and :math:`c_S = c_{0,0}`,
the Monod equation can be expressed in the same form as the Michaelis-Menten kinetics.

Prefactorial extension
^^^^^^^^^^^^^^^^^^^^^^

In some cases a different component, like biomass, may have a prefactorial effect on the reaction rate, which can be represented by a prefactorial extension of the Michaelis-Menten kinetics:

.. math::

    \begin{aligned}
        \nu_{i,j} = \prod_{X = 1}^{N_{X,j}} K_{X,j} c_{X,j} \cdot v_{\mathrm{max},j} \prod_{i = 1}^{N_{sub,j}} \frac{ c_{i,j}}{K_{\mathrm{M}_{i,j}} + c_{i,j}}
    \end{aligned}
   

where:

- :math:`c_{X,j}` is biomass concentration (or any other component that has a prefactorial effect on the reaction rate)
- :math:`K_{X,j}` is the prefactorial constant (if no prefactorial effect it is set as :math:`0`)
- :math:`N_{X,j}` components in :math:`S` that are not in :math:`N_{sub,j}` but have a prefactorial effect on the reaction rate

While those prefactorial components are not contributing to the reaction rate in the same way as substrates, they can still influence the overall reaction rate by acting as any form of inhibitor:

.. math::

    \begin{aligned}
        \nu_{i,j} = \prod_{X = 1}^{N_{X,j}} K_{X,j} c_{X,j} \cdot v_{\mathrm{max},j} \prod_{i = 1}^{N_{sub,j}} \frac{c_{i,j}}{ K_{\mathrm{M}_{i,j}} \,(1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{K^{c}_{I_{k}}}) + c_{i,j} \,(1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{K^{uc}_{I_{k}}})},
    \end{aligned}

where:


- :math:`N_{X,j}` :math:`\cap` :math:`N_{sub,j} = 0`, i.e. the prefactorial components are not substrates
- :math:`N_{X,j}` :math:`\cap` :math:`\mathcal{I}_{i,j} \neq 0`, i.e. the prefactorial components **can** act as inhibitors

Hill Kinetics
^^^^^^^^^^^^^

For substrates without any form of inhibition, CADET supports a Hill-kinetics
extension of the Michaelis-Menten kinetic, introducing a cooperativity exponent :math:`n_{i,j}` for substrate :math:`i` in reaction :math:`j`:

.. math::

    \begin{aligned}
        \nu_{i,j} = \frac{c_{i,j}^{\,n_{i,j}}}{K_{\mathrm{M}_{i,j}}^{\,n_{i,j}} + c_{i,j}^{\,n_{i,j}}}
    \end{aligned}

where

- :math:`n_{i,j}` is the Hill coefficient for substrate :math:`i` in reaction :math:`j`,
- and :math:`K_{\mathrm{M}_{i,j}}`, :math:`c_{i,j}` are defined as above.

For :math:`n_{i,j} = 1`, this reduces to standard Michaelis-Menten kinetics.
Values :math:`n_{i,j} > 1` describe positive cooperativity,
while :math:`0 < n_{i,j} < 1` describes negative cooperativity.

.. note::
    The Hill exponent is only applied to substrates that have **no** competitive or uncompetitive
    inhibitors acting on them (:math:`\mathcal{I}^{c}_{i,j} = \mathcal{I}^{uc}_{i,j} = \emptyset`).

Literature
^^^^^^^^^^

- Segel, I. H. (1993). Enzyme kinetics: Behavior and analysis of rapid equilibrium and steady-state enzyme systems. John Wiley & Sons.
- Monod, Jacques. 1949. “The Growth of Bacterial Cultures.” *Annual Review of Microbiology* 3 (1): 371–394. https://doi.org/10.1146/annurev.mi.03.100149.002103
