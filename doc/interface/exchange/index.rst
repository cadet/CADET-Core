.. _exchange_models:

Exchange Models
===============

Exchange models describe the inter-channel transport processes in the Multi-Channel Transport (MCT) model. These models define how components are exchanged between different channels based on various kinetic and thermodynamic principles.

Available Exchange Models
~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 2

   langmuir_exchange

The exchange models are configured in the MCT unit operation under the ``exchange`` group. Each model implements different kinetic laws for the exchange process:

- **Linear Exchange**: Simple linear exchange kinetics (default)
- **Langmuir Exchange**: Saturation-limited exchange with competitive multi-component kinetics

