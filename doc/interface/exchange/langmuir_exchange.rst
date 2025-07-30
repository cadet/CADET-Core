.. _langmuir_exchange_config:

Langmuir Exchange
=================

**Group /input/model/unit_XXX/exchange – EXCHANGE_MODEL = LANGMUIR_EX**

For information on model equations, refer to :ref:`langmuir_exchange_model`.

The Langmuir Exchange model combines linear exchange kinetics with Langmuir-type capacity limitations for inter-channel transport in the Multi-Channel Transport (MCT) model. 
This model is particularly useful for describing competitive exchange processes where the exchange rate depends on the available capacity in the destination channel.

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

``CAPACITY_MATRIX``
   Maximum exchange capacity matrix

   **Unit:** :math:`mol~m_{channel}^{-3}`

   Matrix containing maximum capacities :math:`q_{max,i}^{comp}` for each component :math:`comp` in each channel :math:`i`. This parameter determines the saturation limit for the Langmuir-type exchange kinetics. The matrix has dimensions :math:`N_{channel} \times N_{comp}`. The vector ordering is channel-component (i-comp) major.

   ===================  =========================  =================================================
   **Type:** double     **Range:** :math:`\gt 0`   **Length:** :math:`N_{channel} \times N_{comp}`
   ===================  =========================  =================================================
