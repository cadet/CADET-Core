.. _exchange_models_theory:

Exchange Models
===============

Exchange models describe the mass transfer kinetics between channels in the Multi-Channel Transport (MCT) model. 
These models define the mathematical relationships that govern how components are exchanged between different phases or compartments.

.. toctree::
   :maxdepth: 2

   langmuir_exchange

Mathematical Framework
~~~~~~~~~~~~~~~~~~~~~~

All exchange models contribute to the MCT mass balance equation through exchange terms:

.. math::

   \frac{\partial c_{i,l}}{\partial t} = \text{Transport Terms} + \sum_{k \neq l} \left(J_{k \rightarrow l}^i - J_{l \rightarrow k}^i\right) + \text{Reaction Terms}

where :math:`J_{l \rightarrow k}^i` represents the exchange flux of component :math:`i` from channel :math:`l` to channel :math:`k`.

The specific form of :math:`J_{l \rightarrow k}^i` depends on the chosen exchange model and can incorporate:

- Concentration driving forces
- Capacity limitations  
- Multi-component competition
