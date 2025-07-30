.. _langmuir_exchange_config:

Langmuir Exchange
=================

**Group /input/model/unit_XXX/exchange – EXCHANGE_MODEL = LANGMUIR_EX**

For information on model equations, refer to :ref:`langmuir_exchange_model`.

The Langmuir Exchange model combines linear exchange kinetics with Langmuir-type saturation limitations for inter-channel transport in the Multi-Channel Transport (MCT) model. 
This model is particularly useful for describing competitive exchange processes where the exchange rate depends on the available saturation level in the destination channel.

``EXCHANGE_MATRIX``
   Exchange rate matrix for inter-channel transport with Langmuir kinetics

   **Unit:** :math:`s^{-1}`

   Matrix containing exchange rate constants :math:`k_{ij}^{comp}` for component :math:`comp` from channel :math:`i` to channel :math:`j`. The matrix has dimensions :math:`N_{channel} \times N_{channel} \times N_{comp}`. The vector ordering is source channel - destination channel - component (i-j-k) major.

   For parameter sensitivity addressing:
   
   - :math:`\texttt{SENS_BOUNDPHASE}` *Channel from* 
   - :math:`\texttt{SENS_PARTYPE}` *Channel to*
   - :math:`\texttt{SENS_COMP}` *Component*

   ===================  =========================  =========================================================================
   **Type:** double     **Range:** :math:`\ge 0`   **Length:** :math:`N_{channel} \times N_{channel} \times N_{comp}`
   ===================  =========================  =========================================================================

``SATURATION_MATRIX``
   Maximum exchange saturation matrix

   **Unit:** :math:`mol~m_{channel}^{-3}`

   Matrix containing maximum saturation levels :math:`q_{max,i}^{comp}` for each component :math:`comp` in each channel :math:`i`. This parameter determines the saturation limit for the Langmuir-type exchange kinetics. The matrix has dimensions :math:`N_{channel} \times N_{comp}`. The vector ordering is channel-component (i-comp) major.

   ===================  =========================  =================================================
   **Type:** double     **Range:** :math:`\gt 0`   **Length:** :math:`N_{channel} \times N_{comp}`
   ===================  =========================  =================================================

``CHANNEL_CROSS_SECTION_AREAS``
   Cross section areas for each channel

   **Unit:** :math:`m^2`

   Defines the cross section area of each channel, which affects the volumetric exchange fluxes between channels. This parameter is inherited from the MCT model configuration but is essential for proper exchange flux calculations.

   ===================  =========================  ==============================
   **Type:** double     **Range:** :math:`\gt 0`   **Length:** :math:`N_{channel}`
   ===================  =========================  ==============================

Model Equations
----------------

The Langmuir Exchange model describes the exchange flux from channel :math:`i` to channel :math:`j` for component :math:`comp` as:

.. math::

   J_{i \rightarrow j}^{comp} = k_{ij}^{comp} \cdot c_i^{comp} \cdot q_{max,j}^{comp} \cdot \left(1 - \sum_{k=1}^{N_{comp}} \frac{c_j^k}{q_{max,j}^k}\right) \cdot \frac{A_i}{A_j}

where:

- :math:`k_{ij}^{comp}` is the exchange rate constant from channel :math:`i` to channel :math:`j` for component :math:`comp`
- :math:`c_i^{comp}` is the concentration of component :math:`comp` in source channel :math:`i`
- :math:`q_{max,j}^{comp}` is the maximum saturation level for component :math:`comp` in destination channel :math:`j`
- :math:`\sum_{k=1}^{N_{comp}} \frac{c_j^k}{q_{max,j}^k}` is the normalized total occupancy in destination channel :math:`j`
- :math:`A_i, A_j` are the cross section areas of channels :math:`i` and :math:`j`

The Langmuir term :math:`\left(1 - \sum_{k=1}^{N_{comp}} \frac{c_j^k}{q_{max,j}^k}\right)` represents the available saturation fraction in the destination channel, which approaches zero as the channel becomes saturated.

Physical Interpretation
-----------------------

The Langmuir Exchange model captures several important physical phenomena:

1. **Competitive Exchange**: The exchange rate decreases as the destination channel becomes more occupied by any component
2. **Saturation Limitation**: Exchange stops when the destination channel reaches its maximum saturation level
3. **Component Competition**: Different components compete for the same exchange sites/saturation levels
4. **Volume Conservation**: Cross section area ratios ensure proper volumetric flux scaling

Applications
------------

This exchange model is particularly suitable for:

- **Liquid-Liquid Chromatography (LLC)**: Modeling competitive partitioning between immiscible phases
- **Multi-Phase Reactors**: Describing mass transfer with saturation constraints
- **Biological Transport**: Modeling saturable transport processes between compartments
- **Adsorption Processes**: Inter-phase transport with competitive adsorption

Example Configuration
---------------------

.. code-block:: json

   {
     "EXCHANGE_MODEL": "LANGMUIR_EX",
     "EXCHANGE_MATRIX": [0.0, 1.5, 2.0, 0.0],
     "SATURATION_MATRIX": [100.0, 150.0, 200.0, 120.0],
     "CHANNEL_CROSS_SECTION_AREAS": [1e-4, 2e-4]
   }

This example configures a 2-channel, 2-component system with:

- Exchange from channel 1 to 2: :math:`k_{12}^1 = 1.5`, :math:`k_{12}^2 = 2.0`
- Exchange from channel 2 to 1: :math:`k_{21}^1 = 0.0`, :math:`k_{21}^2 = 0.0`
- Saturation levels: :math:`q_{max,1}^1 = 100`, :math:`q_{max,1}^2 = 150`, :math:`q_{max,2}^1 = 200`, :math:`q_{max,2}^2 = 120`
- Cross sections: :math:`A_1 = 1 \times 10^{-4}`, :math:`A_2 = 2 \times 10^{-4}` m²

Notes
-----

- The exchange matrix diagonal elements should typically be zero (no self-exchange)
- Saturation values must be positive to avoid numerical instabilities
- The model automatically handles the case where saturation approaches zero by setting exchange to zero
- Numerical stability is ensured by clamping the occupancy sum to the range [0,1]
